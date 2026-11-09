#pragma once

#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include <QGraphicsItem>
#include <QPainterPath>
#include <QPointF>
#include <QRectF>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/Geometry/point.h>

namespace RDKit
{
class Atom;
}

namespace schrodinger::sketcher
{

class AtomDisplaySettings;
class BondDisplaySettings;
class Fonts;

/** Natural 2D size of the atomistic monomer, in molecular coordinates. */
std::optional<RDGeom::Point3D>
get_monomer_atomistic_preview_size(const RDKit::Atom* monomer);

/** Graphics-only atomistic depiction shown while hovering a monomer.
 * Prototype: it does not modify the model or its coordinates.
 */
class MonomerAtomisticPreview final : public QGraphicsItem
{
  public:
    struct AttachmentDirection {
        unsigned int number;
        QPointF position;
        QRectF covered_area;
    };
    using AttachmentDirections = std::vector<AttachmentDirection>;
    MonomerAtomisticPreview(const RDKit::Atom* monomer,
                            const AttachmentDirections& attachment_directions,
                            const Fonts& fonts,
                            const AtomDisplaySettings& atom_settings,
                            const BondDisplaySettings& bond_settings,
                            bool is_dark_mode);
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem*,
               QWidget*) override;

  private:
    std::unique_ptr<RDKit::RWMol> m_mol;
    QRectF m_bounds;
    QPainterPath m_visible_shape;
};

} // namespace schrodinger::sketcher
