#include "schrodinger/sketcher/dialog/sketcher_welcome_dialog.h"
#include "schrodinger/sketcher/ui/ui_sketcher_welcome_dialog.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWelcomeDialog::SketcherWelcomeDialog(QWidget* parent) :
    ModalDialog(parent)
{
    m_ui.reset(new Ui::SketcherWelcomeDialog());
    setupDialogUI(*m_ui);

    connect(m_ui->ok_btn, &QPushButton::clicked, this,
            &SketcherWelcomeDialog::accept);
    connect(m_ui->do_not_show_cb, &QCheckBox::stateChanged, this,
            &SketcherWelcomeDialog::doNotShowCheckboxStateChanged);
}

SketcherWelcomeDialog::~SketcherWelcomeDialog() = default;

void SketcherWelcomeDialog::setDoNotShowCheckboxVisible(bool visible)
{
    m_ui->do_not_show_cb->setVisible(visible);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/sketcher_welcome_dialog.moc"
