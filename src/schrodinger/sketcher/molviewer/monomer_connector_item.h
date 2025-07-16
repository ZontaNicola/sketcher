#pragma once

#include <QGraphicsItem>
#include <QLineF>
#include <QPen>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

class QPainter;

namespace RDKit
{
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class AbstractMonomerItem;

/**
 * A Qt graphics item for representing a connector between two monomers in a
 * molviewer Scene.
 */
class SKETCHER_API MonomerConnectorItem : public AbstractBondOrConnectorItem
{
  public:
    /**
     * Note that this class does not take ownership of bond, so it must not be
     * destroyed while the graphics item is in use.
     *
     * @param bond The RDKit bond that this item should represent.
     *
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     * @pre bond != nullptr
     * @pre bond->hasOwningMol()
     */
    MonomerConnectorItem(const RDKit::Bond* bond,
                         const AbstractMonomerItem& start_monomer_item,
                         const AbstractMonomerItem& end_monomer_item,
                         QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::MONOMER_CONNECTOR) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

  protected:
    QPen m_connector_pen;
    QPen m_arrowhead_pen;
    QBrush m_arrowhead_brush;
    QLineF m_connector_line;
    QPainterPath m_arrowhead_path;
    const AbstractMonomerItem& m_start_item;
    const AbstractMonomerItem& m_end_item;
};

} // namespace sketcher
} // namespace schrodinger
