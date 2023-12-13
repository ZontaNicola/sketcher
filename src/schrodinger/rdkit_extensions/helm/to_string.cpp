#include "schrodinger/rdkit_extensions/helm/to_string.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fmt/format.h>
#include <iterator>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <utility>

#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include "schrodinger/rdkit_extensions/helm.h"

namespace helm
{
namespace
{
[[nodiscard]] std::string
get_polymer_helm(const ::RDKit::SubstanceGroup& polymer);

[[nodiscard]] std::string get_connection_helm(const ::RDKit::Bond* bond);

template <class PropType, class RDKitObject> PropType get_property
    [[nodiscard]] (const RDKitObject* obj, const std::string& propname)
{
    return obj->template getProp<PropType>(propname);
}
} // namespace

[[nodiscard]] std::string rdkit_to_helm(const ::RDKit::ROMol& mol)
{
    if (!schrodinger::rdkit_extensions::is_coarse_grain_mol(mol)) {
        throw std::invalid_argument(
            "Atomistic conversions are currently unsupported");
    }

    const auto& sgroups = ::RDKit::getSubstanceGroups(mol);
    std::vector<::RDKit::SubstanceGroup> polymers;
    std::copy_if(sgroups.begin(), sgroups.end(), std::back_inserter(polymers),
                 [](const auto& sgroup) {
                     return get_property<std::string>(&sgroup, "TYPE") == "COP";
                 });

    std::vector<std::string> polymers_helm;
    std::vector<unsigned int> intra_polymer_bonds;
    std::transform(polymers.begin(), polymers.end(),
                   std::back_inserter(polymers_helm), [&](const auto& polymer) {
                       const auto& bonds = polymer.getBonds();
                       std::copy(bonds.begin(), bonds.end(),
                                 std::back_inserter(intra_polymer_bonds));
                       return get_polymer_helm(polymer);
                   });

    std::stringstream ss;
    ss << boost::join(polymers_helm, "|");

    std::sort(intra_polymer_bonds.begin(), intra_polymer_bonds.end(),
              std::greater<unsigned int>());
    ::RDKit::RWMol mutable_mol(mol);
    for (const auto idx : intra_polymer_bonds) {
        const auto* bond = mutable_mol.getBondWithIdx(idx);
        mutable_mol.removeBond(bond->getBeginAtom()->getIdx(),
                               bond->getEndAtom()->getIdx());
    }
    std::vector<std::string> connections_helm;
    for (size_t idx = 0; idx < mutable_mol.getNumBonds(); ++idx) {
        connections_helm.push_back(
            get_connection_helm(mutable_mol.getBondWithIdx(idx)));
    }

    ss << "$" << boost::join(connections_helm, "|");
    const auto& sgroup =
        *std::find_if(sgroups.begin(), sgroups.end(), [](const auto& sgroup) {
            return get_property<std::string>(&sgroup, "TYPE") == "DAT" &&
                   get_property<std::string>(&sgroup, "FIELDNAME") ==
                       SUPPLEMENTARY_INFORMATION;
        });

    const auto& supplementary_info =
        get_property<std::vector<std::string>>(&sgroup, "DATAFIELDS");
    ss << "$" << supplementary_info[0];
    ss << "$" << supplementary_info[1];
    ss << "$V2.0";
    return ss.str();
}

namespace
{

[[nodiscard]] std::string get_polymer_id(const ::RDKit::Atom* atom)
{
    return static_cast<const RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo())
        ->getChainId();
}

[[nodiscard]] std::string get_monomer_id(const ::RDKit::Atom* atom)
{
    std::string monomer_id;
    if (atom->hasProp(MONOMER_LIST)) {
        monomer_id =
            fmt::format("({})", get_property<std::string>(atom, MONOMER_LIST));
    } else {
        const auto id = get_property<std::string>(atom, ATOM_LABEL);
        const auto& is_blob = (get_polymer_id(atom).front() == 'B');
        monomer_id =
            ((id.size() == 1 || is_blob) ? id : fmt::format("[{}]", id));
    }
    if (atom->hasProp(ANNOTATION)) {
        monomer_id.append('"' + get_property<std::string>(atom, ANNOTATION) +
                          '"');
    }
    return monomer_id;
}

[[nodiscard]] std::string
get_polymer_helm(const ::RDKit::SubstanceGroup& polymer)
{
    const auto& mol = polymer.getOwningMol();
    std::unordered_set<int> branch_monomers;
    const auto bonds = polymer.getBonds();
    std::for_each(bonds.begin(), bonds.end(), [&](const auto& idx) {
        const auto bond = mol.getBondWithIdx(idx);
        const auto atom2 = bond->getEndAtom();
        if (!atom2->hasProp(ATOM_LABEL)) { // for the temp query atom
            return;
        }
        if (get_property<std::string>(bond, LINKAGE) == BRANCH_LINKAGE) {
            branch_monomers.insert(atom2->getIdx());
        }
    });

    std::unordered_map<int, int> repeated_group_starts;
    const auto& sgroups = ::RDKit::getSubstanceGroups(mol);
    std::for_each(sgroups.begin(), sgroups.end(), [&](const auto& sgroup) {
        if (get_property<std::string>(&sgroup, "TYPE") != "SRU") {
            return;
        }
        repeated_group_starts.insert(
            {sgroup.getAtoms().front(), sgroup.getIndexInMol()});
    });

    std::stringstream ss;

    int sru_idx = -1;
    bool prev_is_branch = false;
    const auto& atoms = polymer.getAtoms();
    ss << get_polymer_id(mol.getAtomWithIdx(atoms.front())) << "{";

    std::for_each(atoms.begin(), atoms.end(), [&](const auto& idx) {
        const auto atom = mol.getAtomWithIdx(idx);
        if (!atom->hasProp(ATOM_LABEL)) { // for the temp query atom
            return;
        }

        if (repeated_group_starts.count(idx)) {
            sru_idx = repeated_group_starts.at(idx);
            ss << (idx == atoms.front() ? "" : ".");
            ss << (sgroups[sru_idx].getAtoms().size() == 1 ? "" : "(");
        }

        if (idx == atoms.front() || prev_is_branch ||
            repeated_group_starts.count(idx)) {
            ss << get_monomer_id(atom);
            prev_is_branch = false;
        } else if (branch_monomers.count(idx)) {
            ss << "(" << get_monomer_id(atom) << ")";
            prev_is_branch = true;
        } else {
            ss << "." << get_monomer_id(atom);
        }

        if (sru_idx != -1 && sgroups[sru_idx].getAtoms().back() == idx) {
            const auto& sgroup = sgroups[sru_idx];
            ss << (sgroup.getAtoms().size() == 1 ? "" : ")");
            ss << "'" << get_property<std::string>(&sgroup, "LABEL") << "'";
            if (sgroup.hasProp(ANNOTATION)) {
                ss << '"' << get_property<std::string>(&sgroup, ANNOTATION)
                   << '"';
            }
            sru_idx = -1;
        }
    });
    ss << "}";
    if (polymer.hasProp(ANNOTATION)) {
        ss << '"' << get_property<std::string>(&polymer, ANNOTATION) << '"';
    }
    return ss.str();
};

[[nodiscard]] std::string get_connection_helm(const ::RDKit::Bond* bond)
{
    // helper function to always make BLOB linkage information ambiguous
    static auto process_linkage_info = [](const auto& polymer_id,
                                          const auto& linkage_info) {
        return polymer_id.front() == 'B' ? "?" : linkage_info;
    };

    static auto get_linkage_residue = [](const auto* atom) {
        auto* res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        auto polymer_id = res_info->getChainId();
        auto resnum = std::to_string(res_info->getResidueNumber());
        return std::make_pair<std::string, std::string>(
            std::move(polymer_id), process_linkage_info(polymer_id, resnum));
    };

    const auto [from_rgroup, to_rgroup] = [&]() {
        auto linkage = get_property<std::string>(bond, LINKAGE);
        // NOTE: linkages can be of the forms RX-RX, RX-?, ?-RX, ?-?,pair-pair
        // etc
        auto middle = linkage.find("-");
        return std::make_pair<std::string, std::string>(
            linkage.substr(0, middle), linkage.substr(middle + 1));
    }();

    auto [from_id, from_res] = get_linkage_residue(bond->getBeginAtom());
    auto [to_id, to_res] = get_linkage_residue(bond->getEndAtom());

    auto connection =
        fmt::format("{},{},{}:{}-{}:{}", from_id, to_id, from_res,
                    process_linkage_info(from_id, from_rgroup), to_res,
                    process_linkage_info(to_id, to_rgroup));

    if (bond->hasProp("ANNOTATION")) {
        connection.append(fmt::format(
            "\"{}\"", get_property<std::string>(bond, "ANNOTATION")));
    }
    return connection;
}
} // namespace
} // namespace helm
