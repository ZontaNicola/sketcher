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

using schrodinger::rdkit_extensions::compute_monomer_mol_coords;
using schrodinger::rdkit_extensions::coordinates_are_valid;
using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::to_rdkit;
using schrodinger::sketcher::get_image_bytes;
using schrodinger::sketcher::ImageFormat;
using schrodinger::sketcher::RenderOptions;

namespace
{

void write_response(const bool success, const bool coordinates_valid,
                    const QByteArray& payload)
{
    // Base64 keeps SVG, HELM, and exception text from introducing tabs or
    // newlines into the one-request-per-line protocol.
    std::cout << success << '\t' << coordinates_valid << '\t'
              << payload.toBase64().toStdString() << '\n';
    std::cout.flush();
}

} // namespace

int main(int argc, char* argv[])
{
    QApplication application(argc, argv);
    RenderOptions options;
    options.width_height = {400, 400};
    options.background_color = Qt::white;

    std::string encoded_helm;
    while (std::getline(std::cin, encoded_helm)) {
        try {
            const auto helm_bytes =
                QByteArray::fromBase64(QByteArray::fromStdString(encoded_helm));
            auto mol = to_rdkit(helm_bytes.toStdString(), Format::HELM);
            compute_monomer_mol_coords(*mol);
            const auto coordinates_valid = coordinates_are_valid(*mol);
            const auto svg = get_image_bytes(*mol, ImageFormat::SVG, options);
            write_response(true, coordinates_valid, svg);
        } catch (const std::exception& error) {
            write_response(false, false, QByteArray(error.what()));
        }
    }
    return 0;
}
