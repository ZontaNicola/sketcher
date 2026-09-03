
#define BOOST_TEST_MODULE monomeric

#include <unordered_set>

#include <rdkit/GraphMol/RWMol.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <QGraphicsRectItem>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

using namespace boost::unit_test;

namespace schrodinger
{
namespace sketcher
{

/**
 * Make sure that contains_two_monomer_linkages correctly detects two monomer
 * linkages in the same bond when there's a disulfide bond between neighboring
 * cysteines.
 */
BOOST_AUTO_TEST_CASE(test_contains_two_monomer_linkages)
{
    // two neighboring cysteines, but no disulfide
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{C.C}$$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(0)));

    // two neighboring cysteines with a disulfide
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(contains_two_monomer_linkages(mol->getBondWithIdx(0)));

    // a disulfide, but between two non-neighboring cysteines
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 3);
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(0)));
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(1)));
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(2)));
}

/**
 * Make sure that get_bound_attachment_point_names_and_atoms() and
 * get_available_attachment_point_names() return the expected attachment point
 * names for a variety of molecules
 */
BOOST_AUTO_TEST_CASE(test_get_attachment_points)
{
    std::vector<BoundAttachmentPoint> bound_aps, exp_bound;
    std::vector<UnboundAttachmentPoint> unbound_aps, exp_available;
    const RDKit::Atom *atom0, *atom1, *atom2;

    // a lone alanine has no bound attachment points
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"R1", "N", 1, Direction::W},
                         {"R2", "C", 2, Direction::E},
                         {"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // two alanines next to each other
    mol = rdkit_extensions::to_rdkit("PEPTIDE1{A.A}$$$$V2.0");
    prepare_mol(*mol);
    {
        const auto* atom0 = mol->getAtomWithIdx(0);
        const auto* atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W},
                         {"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E},
                         {"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // a side-chain interaction between two non-adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0");
    prepare_mol(*mol);
    {
        // put the three residues in a horizontal line
        auto& conf = mol->getConformer();
        conf.setAtomPos(0, {-BOND_LENGTH, 0.0, 0.0});
        conf.setAtomPos(1, {0.0, 0.0, 0.0});
        conf.setAtomPos(2, {BOND_LENGTH, 0.0, 0.0});

        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);
        atom2 = mol->getAtomWithIdx(2);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R3", "S", 3, atom2, false, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R2", "C", 2, atom2, false, Direction::E}};
        exp_available = {{"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom2);
        exp_bound = {{"R1", "N", 1, atom1, false, Direction::W},
                     {"R3", "S", 3, atom0, false, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // a side-chain interaction between two adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R3", "S", 3, atom1, true, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R3", "S", 3, atom0, true, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // CHEM monomers
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[MONO1]}|CHEM2{[MONO2]}$CHEM1,CHEM2,1:R1-1:R3$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", "R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", "R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R3", "R3", 3, atom0, false, Direction::N}};
        exp_available = {{"R1", "R1", 1, Direction::W},
                         {"R2", "R2", 2, Direction::E},
                         {"R4", "R4", 4, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // CHEM monomers with too many attachment points - some will be omitted
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[MONO1]}|CHEM2{[MONO2]}$CHEM1,CHEM2,1:R1-1:R11$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", "R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", "R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R11", "R11", 11, atom0, false, Direction::N}};
        exp_available = {
            {"R1", "R1", 1, Direction::W},  {"R2", "R2", 2, Direction::E},
            {"R3", "R3", 3, Direction::S},  {"R4", "R4", 4, Direction::NW},
            {"R5", "R5", 5, Direction::NE}, {"R6", "R6", 6, Direction::SE},
            {"R7", "R7", 7, Direction::SW}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // NA_PHOSPHATE monomers name their unbound attachment points after the
    // attachment point of the bound sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P}$$$$");
    prepare_mol(*mol);
    {
        auto start_phos = mol->getAtomWithIdx(0);
        auto start_sugar = mol->getAtomWithIdx(1);
        auto middle_phos = mol->getAtomWithIdx(3);
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phosphate = mol->getAtomWithIdx(6);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(start_phos);
        exp_bound = {{"R2", "", 2, start_sugar, false, Direction::E}};
        exp_available = {{"R1", "5'", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(middle_phos);
        exp_bound = {{"R1", "", 1, start_sugar, false, Direction::W},
                     {"R2", "", 2, term_sugar, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phosphate);
        exp_bound = {{"R1", "", 1, term_sugar, false, Direction::W}};
        exp_available = {{"R2", "3'", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // NA_PHOSPHATE monomers still take their name from the bound sugar even if
    // they're at the end of a chain of phosphates
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P.P.P}$$$$");
    prepare_mol(*mol);
    {
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phos_chain_1 = mol->getAtomWithIdx(6);
        auto term_phos_chain_2 = mol->getAtomWithIdx(7);
        auto term_phos_chain_3 = mol->getAtomWithIdx(8);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_1);
        exp_bound = {{"R1", "", 1, term_sugar, false, Direction::W},
                     {"R2", "", 2, term_phos_chain_2, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_2);
        exp_bound = {{"R1", "", 1, term_phos_chain_1, false, Direction::W},
                     {"R2", "", 2, term_phos_chain_3, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_3);
        exp_bound = {{"R1", "", 1, term_phos_chain_2, false, Direction::W}};
        exp_available = {{"R2", "3'", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // attachment points on lone phosphates are unnamed since there's no bound
    // sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P}$$$$");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"R1", "", 1, Direction::W},
                         {"R2", "", 2, Direction::E}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // an amino acid with unrecognized attachment points
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{A.A}$PEPTIDE1,PEPTIDE1,1:R4-2:R4$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R4", "R4", 4, atom1, true, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W},
                         {"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R4", "R4", 4, atom0, true, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E},
                         {"pair", "H-bond", -1, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // single stranded RNA
    mol = rdkit_extensions::to_rdkit("RNA1{R(A)P.R(C)P.R(G)P}$$$$V2.0");
    prepare_mol(*mol);
    {
        auto sugar = mol->getAtomWithIdx(0);
        auto base = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(base);
        exp_bound = {{"R1", "N1/9", 1, sugar, false, Direction::N}};
        exp_available = {{"pair", "H-bond", ATTACHMENT_POINT_WITH_CUSTOM_NAME,
                          Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // double stranded RNA
    mol = rdkit_extensions::to_rdkit(
        "RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}"
        "$RNA1,RNA2,2:pair-9:pair"
        "|RNA1,RNA2,5:pair-6:pair"
        "|RNA1,RNA2,8:pair-3:pair$$$V2.0");
    prepare_mol(*mol);
    {
        auto sugar = mol->getAtomWithIdx(0);
        auto base = mol->getAtomWithIdx(1);
        auto paired_base = mol->getAtomWithIdx(17);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(base);
        exp_bound = {{"R1", "N1/9", 1, sugar, false, Direction::N},
                     {"pair", "pair", ATTACHMENT_POINT_WITH_CUSTOM_NAME,
                      paired_base, false, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());
    }
}

/**
 * Make sure that get_attachment_points returns the correct number for parseable
 * chain name, and returns -1 when the chain name can't be parsed.
 */
BOOST_AUTO_TEST_CASE(test_get_chain_num)
{
    using rdkit_extensions::ChainType;
    BOOST_TEST(get_chain_num("PEPTIDE1", ChainType::PEPTIDE) == 1);
    BOOST_TEST(get_chain_num("PEPTIDE3", ChainType::PEPTIDE) == 3);
    BOOST_TEST(get_chain_num("PEPTIDE", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("ABCDEFG", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("ABCDEFG2", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("PEPTIDE1", ChainType::RNA) == -1);
    BOOST_TEST(get_chain_num("RNA2", ChainType::RNA) == 2);
}

/**
 * Make sure that get_first_available_chain_name returns the expected chain name
 * of the appropriate chain type.
 */
BOOST_AUTO_TEST_CASE(test_get_first_available_chain_name)
{
    using rdkit_extensions::ChainType;

    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}|PEPTIDE2{A}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE3");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA1");

    mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}|PEPTIDE3{A}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE2");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA1");

    mol = rdkit_extensions::to_rdkit(
        "RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA3");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE1");
}

BOOST_AUTO_TEST_CASE(test_monomer_arrowhead_offset_uses_available_side)
{
    QGraphicsRectItem monomer(-10.0, -5.0, 20.0, 10.0);
    monomer.setPos(100.0, 100.0);

    const auto right_offset =
        get_monomer_arrowhead_offset(monomer, QPointF(200.0, 100.0), {});
    BOOST_TEST(right_offset.x() == monomer.boundingRect().right() +
                                       MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(right_offset.y() == 0.0, boost::test_tools::tolerance(0.001));

    const auto left_offset =
        get_monomer_arrowhead_offset(monomer, QPointF(0.0, 100.0), {});
    BOOST_TEST(left_offset.x() == monomer.boundingRect().left() -
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(left_offset.y() == 0.0, boost::test_tools::tolerance(0.001));

    const auto bottom_offset =
        get_monomer_arrowhead_offset(monomer, QPointF(100.0, 200.0), {});
    BOOST_TEST(bottom_offset.x() == 0.0, boost::test_tools::tolerance(0.001));
    BOOST_TEST(bottom_offset.y() == monomer.boundingRect().bottom() +
                                        MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));

    const auto fallback_offset = get_monomer_arrowhead_offset(
        monomer, QPointF(200.0, 50.0), {Direction::E});
    BOOST_TEST(fallback_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(fallback_offset.y() == monomer.boundingRect().top() -
                                          MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));

    const auto corner_offset = get_monomer_arrowhead_offset(
        monomer, QPointF(200.0, 50.0), {Direction::E, Direction::N});
    QLineF expected_corner_offset(QPointF(), monomer.boundingRect().topRight());
    expected_corner_offset.setLength(expected_corner_offset.length() +
                                     MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    BOOST_TEST(corner_offset.x() == expected_corner_offset.p2().x(),
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(corner_offset.y() == expected_corner_offset.p2().y(),
               boost::test_tools::tolerance(0.001));
}

BOOST_AUTO_TEST_CASE(test_duplicate_custom_bond_does_not_occupy_side)
{
    auto mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{[dC].[dQ].I.[dD].S.[dP].[dC]}$PEPTIDE1,PEPTIDE1,"
        "1:R3-7:R3$$$");
    prepare_mol(*mol);
    auto& conf = mol->getConformer();
    conf.setAtomPos(0, {0.0, 0.0, 0.0});
    conf.setAtomPos(1, {BOND_LENGTH, 0.0, 0.0});
    conf.setAtomPos(5, {BOND_LENGTH, -BOND_LENGTH, 0.0});
    conf.setAtomPos(6, {0.0, -BOND_LENGTH, 0.0});

    const auto* first_monomer = mol->getAtomWithIdx(0);
    const auto* last_monomer = mol->getAtomWithIdx(6);
    QGraphicsRectItem first_item(-10.0, -5.0, 20.0, 10.0);
    QGraphicsRectItem last_item(-10.0, -5.0, 20.0, 10.0);
    last_item.setPos(0.0, BOND_LENGTH);

    const auto first_offset = get_monomer_arrowhead_offset(
        first_item, last_item.pos(), first_monomer, last_monomer, false);
    BOOST_TEST(first_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(first_offset.y() == first_item.boundingRect().bottom() +
                                       MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));

    const auto last_offset = get_monomer_arrowhead_offset(
        last_item, first_item.pos(), last_monomer, first_monomer, false);
    BOOST_TEST(last_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(last_offset.y() == last_item.boundingRect().top() -
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));
}

BOOST_AUTO_TEST_CASE(test_diagonal_connections_do_not_occupy_sides)
{
    auto mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.K.G.K.G.A.K.C.S.R.L.M.Y.D.C.C.T.G.S.C.R.S.G.K.C}"
        "$PEPTIDE1,PEPTIDE1,1:R3-16:R3|PEPTIDE1,PEPTIDE1,8:R3-20:R3|"
        "PEPTIDE1,PEPTIDE1,15:R3-25:R3$$$");
    rdkit_extensions::compute_monomer_mol_coords(*mol);

    const auto* monomer_8 = mol->getAtomWithIdx(7);
    const auto* monomer_20 = mol->getAtomWithIdx(19);
    const auto& conf = mol->getConformer();
    auto monomer_8_pos = to_scene_xy(conf.getAtomPos(monomer_8->getIdx()));
    auto monomer_20_pos = to_scene_xy(conf.getAtomPos(monomer_20->getIdx()));

    QGraphicsRectItem monomer_8_item(-10.0, -5.0, 20.0, 10.0);
    monomer_8_item.setPos(monomer_8_pos);
    QGraphicsRectItem monomer_20_item(-10.0, -5.0, 20.0, 10.0);
    monomer_20_item.setPos(monomer_20_pos);

    const auto monomer_8_offset = get_monomer_arrowhead_offset(
        monomer_8_item, monomer_20_pos, monomer_8, monomer_20, false);
    BOOST_TEST(monomer_8_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(monomer_8_offset.y() ==
                   monomer_8_item.boundingRect().top() -
                       MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));

    const auto monomer_20_offset = get_monomer_arrowhead_offset(
        monomer_20_item, monomer_8_pos, monomer_20, monomer_8, false);
    BOOST_TEST(monomer_20_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(monomer_20_offset.y() ==
                   monomer_20_item.boundingRect().bottom() +
                       MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
               boost::test_tools::tolerance(0.001));
}

BOOST_AUTO_TEST_CASE(test_parallel_connections_use_same_fallback_side)
{
    auto mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$");
    prepare_mol(*mol);
    auto& conf = mol->getConformer();
    conf.setAtomPos(0, {0.0, 0.0, 0.0});
    // A small vertical perturbation used to send the two arrowheads to
    // opposite sides, causing the disulfide connector to cross the backbone.
    conf.setAtomPos(1, {BOND_LENGTH, 0.1 * BOND_LENGTH, 0.0});

    const auto* first_monomer = mol->getAtomWithIdx(0);
    const auto* second_monomer = mol->getAtomWithIdx(1);
    auto first_pos = to_scene_xy(conf.getAtomPos(0));
    auto second_pos = to_scene_xy(conf.getAtomPos(1));
    QGraphicsRectItem first_item(-10.0, -5.0, 20.0, 10.0);
    first_item.setPos(first_pos);
    QGraphicsRectItem second_item(-10.0, -5.0, 20.0, 10.0);
    second_item.setPos(second_pos);

    const auto first_offset = get_monomer_arrowhead_offset(
        first_item, second_pos, first_monomer, second_monomer, true);
    const auto second_offset = get_monomer_arrowhead_offset(
        second_item, first_pos, second_monomer, first_monomer, true);
    const auto expected_y = first_item.boundingRect().top() -
                            MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
    BOOST_TEST(first_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(second_offset.x() == 0.0,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(first_offset.y() == expected_y,
               boost::test_tools::tolerance(0.001));
    BOOST_TEST(second_offset.y() == expected_y,
               boost::test_tools::tolerance(0.001));
}

} // namespace sketcher
} // namespace schrodinger
