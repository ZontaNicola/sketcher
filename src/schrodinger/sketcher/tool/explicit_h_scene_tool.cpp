#include "schrodinger/sketcher/tool/explicit_h_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"

namespace schrodinger
{
namespace sketcher
{

ExplicitHsSceneTool::ExplicitHsSceneTool(Scene* scene, MolModel* mol_model) :
    SceneToolWithPredictiveHighlighting(scene, mol_model)
{
    m_highlight_types = InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
}

void ExplicitHsSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        auto atom = atom_item->getAtom();
        auto action =
            (m_mol_model->hasAnyImplicitHs({atom}) ? ExplicitHActions::ADD
                                                   : ExplicitHActions::REMOVE);
        m_mol_model->updateExplicitHs(action, {atom});
    }
}

} // namespace sketcher
} // namespace schrodinger
