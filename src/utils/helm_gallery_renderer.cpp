/* -------------------------------------------------------------------------
 * Render HELM strings to SVG using the locally built Sketcher libraries.
 *
 * Reads one base64-encoded HELM string per line. Writes tab-separated status,
 * coordinate-validity, and base64-encoded SVG/error payload fields.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include <iostream>
#include <string>

#include <QApplication>
#include <QByteArray>
#include <QColor>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/sketcher/image_generation.h"

using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::compute_monomer_mol_coords;
using schrodinger::rdkit_extensions::coordinates_are_valid;
using schrodinger::rdkit_extensions::to_rdkit;
using schrodinger::sketcher::ImageFormat;
using schrodinger::sketcher::RenderOptions;
using schrodinger::sketcher::get_image_bytes;

int main(int argc, char* argv[])
{
    QApplication application(argc, argv);
    RenderOptions options;
    options.width_height = {400, 400};
    options.background_color = Qt::white;

    std::string encoded_helm;
    while (std::getline(std::cin, encoded_helm)) {
        try {
            auto helm_bytes = QByteArray::fromBase64(
                QByteArray::fromStdString(encoded_helm));
            auto helm = helm_bytes.toStdString();
            auto mol = to_rdkit(helm, Format::HELM);
            compute_monomer_mol_coords(*mol);
            auto coords_valid = coordinates_are_valid(*mol);
            auto svg = get_image_bytes(*mol, ImageFormat::SVG, options);
            std::cout << "1\t" << coords_valid << "\t"
                      << svg.toBase64().toStdString() << '\n';
        } catch (const std::exception& error) {
            auto encoded_error =
                QByteArray(error.what()).toBase64().toStdString();
            std::cout << "0\t0\t" << encoded_error << '\n';
        }
        std::cout.flush();
    }
    return 0;
}
