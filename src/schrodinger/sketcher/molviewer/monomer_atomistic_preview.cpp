#include "schrodinger/sketcher/molviewer/monomer_atomistic_preview.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <unordered_map>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger::sketcher
{

namespace
{

constexpr double PREVIEW_SCALE = 0.66;

std::unique_ptr<RDKit::RWMol> make_preview_mol(
    const RDKit::Atom* monomer,
    std::unordered_map<unsigned int, unsigned int>* mapped_atoms = nullptr)
{
    const auto label = monomer->getProp<std::string>(ATOM_LABEL);
    bool inline_smiles = false;
    monomer->getPropIfPresent(SMILES_MONOMER, inline_smiles);
    const std::optional<std::string> smiles =
        inline_smiles
            ? std::optional<std::string>(label)
            : rdkit_extensions::MonomerDatabase::instance().getMonomerSmiles(
                  label, rdkit_extensions::getChainType(*monomer));
    if (!smiles)
        return nullptr;

    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(*smiles, 0, false));
    if (!mol || mol->getNumAtoms() == 0)
        return nullptr;
    for (auto* atom : mol->atoms()) {
        std::string rgroup_label;
        if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                   rgroup_label) &&
            rgroup_label.starts_with("_R")) {
            atom->setAtomMapNum(std::stoi(rgroup_label.substr(2)));
            atom->clearProp(RDKit::common_properties::atomLabel);
        }
        if (atom->getAtomMapNum() != 0) {
            if (mapped_atoms != nullptr)
                (*mapped_atoms)[atom->getAtomMapNum()] = atom->getIdx();
            // A carbon stands in for the atom of the connected residue.
            atom->setAtomicNum(6);
            atom->setAtomMapNum(0);
            atom->setFormalCharge(0);
            atom->setNumExplicitHs(0);
            atom->setNoImplicit(false);
        }
    }
    rdkit_extensions::compute2DCoords(*mol);
    return mol;
}

} // namespace

std::optional<RDGeom::Point3D>
get_monomer_atomistic_preview_size(const RDKit::Atom* monomer)
{
    try {
        auto mol = make_preview_mol(monomer);
        if (!mol)
            return std::nullopt;
        const auto& conf = mol->getConformer();
        double min_x = conf.getAtomPos(0).x, max_x = min_x;
        double min_y = conf.getAtomPos(0).y, max_y = min_y;
        for (const auto* atom : mol->atoms()) {
            const auto p = conf.getAtomPos(atom->getIdx());
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
        }
        return RDGeom::Point3D((max_x - min_x + BOND_LENGTH) * PREVIEW_SCALE,
                               (max_y - min_y + BOND_LENGTH) * PREVIEW_SCALE,
                               0);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

MonomerAtomisticPreview::MonomerAtomisticPreview(
    const RDKit::Atom* monomer, const AttachmentDirections& directions,
    const Fonts& fonts, const AtomDisplaySettings& atom_settings,
    const BondDisplaySettings& bond_settings, bool is_dark_mode)
{
    std::unordered_map<unsigned int, unsigned int> mapped_atoms;
    m_mol = make_preview_mol(monomer, &mapped_atoms);
    if (!m_mol)
        return;
    auto& conf = m_mol->getConformer();
    double min_x = conf.getAtomPos(0).x, max_x = min_x;
    double min_y = conf.getAtomPos(0).y, max_y = min_y;
    for (const auto* atom : m_mol->atoms()) {
        const auto p = conf.getAtomPos(atom->getIdx());
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }
    const QPointF center((min_x + max_x) / 2, -(min_y + max_y) / 2);
    auto scene_point = [&](unsigned int index) {
        const auto p = conf.getAtomPos(index);
        return (QPointF(p.x, -p.y) - center) * (VIEW_SCALE * PREVIEW_SCALE);
    };

    double cos_sum = 0, sin_sum = 0;
    unsigned int matches = 0;
    for (const auto& direction : directions) {
        const auto found = mapped_atoms.find(direction.number);
        if (found == mapped_atoms.end())
            continue;
        const auto source = scene_point(found->second);
        if (qFuzzyIsNull(QPointF::dotProduct(source, source)))
            continue;
        const auto angle =
            std::atan2(direction.position.y(), direction.position.x()) -
            std::atan2(source.y(), source.x());
        cos_sum += std::cos(angle);
        sin_sum += std::sin(angle);
        ++matches;
    }
    const auto angle =
        matches && (!qFuzzyIsNull(cos_sum) || !qFuzzyIsNull(sin_sum))
            ? std::atan2(sin_sum, cos_sum)
            : 0.0;
    const auto cosine = std::cos(angle), sine = std::sin(angle);
    for (const auto* atom : m_mol->atoms()) {
        const auto p = scene_point(atom->getIdx());
        const QPointF transformed(p.x() * cosine - p.y() * sine,
                                  p.x() * sine + p.y() * cosine);
        conf.setAtomPos(atom->getIdx(), to_mol_xy(transformed));
    }
    // Anchor each temporary carbon exactly at its connected residue. The
    // residue's graphics will cover it and the last part of the bond.
    for (const auto& direction : directions) {
        const auto found = mapped_atoms.find(direction.number);
        if (found != mapped_atoms.end())
            conf.setAtomPos(found->second, to_mol_xy(direction.position));
    }

    auto [items, atoms, bonds, secondary, sgroups] =
        create_graphics_items_for_mol(m_mol.get(), fonts, atom_settings,
                                      bond_settings, is_dark_mode);
    for (auto* item : items)
        item->setParentItem(this);
    m_bounds = childrenBoundingRect();
    m_visible_shape.addRect(m_bounds);
    for (const auto& direction : directions) {
        QPainterPath covered;
        covered.addRect(direction.covered_area);
        m_visible_shape = m_visible_shape.subtracted(covered);
    }
    setFlag(QGraphicsItem::ItemClipsChildrenToShape);
    setZValue(static_cast<qreal>(ZOrder::MONOMER) + 0.5);
    setAcceptedMouseButtons(Qt::NoButton);
}

QRectF MonomerAtomisticPreview::boundingRect() const
{
    return m_bounds;
}
QPainterPath MonomerAtomisticPreview::shape() const
{
    return m_visible_shape;
}
void MonomerAtomisticPreview::paint(QPainter*, const QStyleOptionGraphicsItem*,
                                    QWidget*)
{
}

} // namespace schrodinger::sketcher
