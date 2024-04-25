#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fasta_to_rdkit

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <string>

#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/fasta_examples.h"
#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"

namespace bdata = boost::unit_test::data;

BOOST_DATA_TEST_CASE(TestInvalidPeptides,
                     bdata::make(fasta::INVALID_PEPTIDE_EXAMPLES), input_fasta)
{
    BOOST_CHECK_THROW(fasta::peptide_fasta_to_rdkit(input_fasta),
                      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(TestInvalidNucleotides,
                     bdata::make(fasta::INVALID_NUCLEOTIDE_EXAMPLES),
                     input_fasta)
{
    BOOST_CHECK_THROW(fasta::rna_fasta_to_rdkit(input_fasta),
                      std::invalid_argument);

    BOOST_CHECK_THROW(fasta::dna_fasta_to_rdkit(input_fasta),
                      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(TestValidPeptides,
                     bdata::make(fasta::VALID_PEPTIDE_EXAMPLES), input_fasta)
{
    BOOST_CHECK_NO_THROW(fasta::peptide_fasta_to_rdkit(input_fasta));
}

BOOST_DATA_TEST_CASE(TestValidNucleotides,
                     bdata::make(fasta::VALID_NUCLEOTIDE_EXAMPLES), input_fasta)
{
    BOOST_CHECK_NO_THROW(fasta::rna_fasta_to_rdkit(input_fasta));
    BOOST_CHECK_NO_THROW(fasta::dna_fasta_to_rdkit(input_fasta));
}

BOOST_DATA_TEST_CASE(
    TestHelmConversionOfPeptideFasta,
    bdata::make(std::vector<std::string>{
        ">\nAAZL", ">\n-A-APL", ">\nAAPL", ">\nARGXCKXEDA",
        ">some description\nAAA", ">some description\nAAA\nDDD",
        ">some description\nAAA\n>\nDDD", ">some description\nAAA\n\nDDD"}) ^
        bdata::make(std::vector<std::string>{
            "PEPTIDE1{A.A.(E+Q).L}$$$$V2.0", "PEPTIDE1{A.A.P.L}$$$$V2.0",
            "PEPTIDE1{A.A.P.L}$$$$V2.0",
            "PEPTIDE1{A.R.G.X.C.K.X.E.D.A}$$$$V2.0",
            R"(PEPTIDE1{A.A.A}"some description"$$$$V2.0)",
            R"(PEPTIDE1{A.A.A.D.D.D}"some description"$$$$V2.0)",
            R"(PEPTIDE1{A.A.A}"some description"|PEPTIDE2{D.D.D}$$$$V2.0)",
            R"(PEPTIDE1{A.A.A.D.D.D}"some description"$$$$V2.0)"}),
    peptide_fasta, equivalent_helm)
{
    auto mol = fasta::peptide_fasta_to_rdkit(peptide_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == equivalent_helm);
}

BOOST_DATA_TEST_CASE(TestHelmConversionOfRNAFasta,
                     bdata::make(std::vector<std::string>{
                         ">\nAAK",
                         ">\nAAB",
                     }) ^ bdata::make(std::vector<std::string>{
                              "RNA1{R(A)P.R(A)P.R((G+T+U))P}$$$$V2.0",
                              "RNA1{R(A)P.R(A)P.R((C+G+T+U))P}$$$$V2.0",
                          }),
                     nucleotide_fasta, equivalent_helm)
{
    auto mol = fasta::rna_fasta_to_rdkit(nucleotide_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == equivalent_helm);
}

BOOST_DATA_TEST_CASE(TestHelmConversionOfDNAFasta,
                     bdata::make(std::vector<std::string>{
                         ">\nAAK",
                         ">\nAAB",
                     }) ^
                         bdata::make(std::vector<std::string>{
                             "RNA1{[dR](A)P.[dR](A)P.[dR]((G+T+U))P}$$$$V2.0",
                             "RNA1{[dR](A)P.[dR](A)P.[dR]((C+G+T+U))P}$$$$V2.0",
                         }),
                     nucleotide_fasta, equivalent_helm)
{
    auto mol = fasta::dna_fasta_to_rdkit(nucleotide_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == equivalent_helm);
}
