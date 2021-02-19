
/*LICENSE_START*/
/*
 *  Copyright (C) 2014 Washington University School of Medicine
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*LICENSE_END*/

#include <QAction>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QLabel>
#include <QStackedWidget>
#include <QToolButton>

#define __BRAIN_BROWSER_WINDOW_TOOL_BAR_CHART_TWO_ATTRIBUTES_DECLARE__
#include "BrainBrowserWindowToolBarChartTwoAttributes.h"
#undef __BRAIN_BROWSER_WINDOW_TOOL_BAR_CHART_TWO_ATTRIBUTES_DECLARE__

#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "CaretDataFile.h"
#include "CaretDataFileSelectionModel.h"
#include "EnumComboBoxTemplate.h"
#include "CaretMappableDataFile.h"
#include "ChartTwoMatrixDisplayProperties.h"
#include "EventGraphicsUpdateAllWindows.h"
#include "EventManager.h"
#include "EventUserInterfaceUpdate.h"
#include "ModelChartTwo.h"
#include "WuQFactory.h"
#include "WuQMacroManager.h"
#include "WuQWidgetObjectGroup.h"
#include "WuQtUtilities.h"

using namespace caret;


    
/**
 * \class caret::BrainBrowserWindowToolBarChartTwoAttributes 
 * \brief Controls for chart attributes.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param parentToolBar
 *   The parent toolbar.
 * @param parentObjectName
 *   Name of parent for macros
 */
BrainBrowserWindowToolBarChartTwoAttributes::BrainBrowserWindowToolBarChartTwoAttributes(BrainBrowserWindowToolBar* parentToolBar,
                                                                                         const QString& parentObjectName)
: BrainBrowserWindowToolBarComponent(parentToolBar)
{
    const QString objectNamePrefix(parentObjectName
                                   + ":ChartTwoAttributes");
    
    m_cartesianChartAttributesWidget = new CartesianChartTwoAttributesWidget(this,
                                                                             objectNamePrefix);
    
    m_matrixChartTwoAttributesWidget = new MatrixChartTwoAttributesWidget(this,
                                                                          objectNamePrefix);
    
    m_stackedWidget = new QStackedWidget();
    m_stackedWidget->addWidget(m_cartesianChartAttributesWidget);
    m_stackedWidget->addWidget(m_matrixChartTwoAttributesWidget);
    
    QVBoxLayout* layout = new QVBoxLayout(this);
    WuQtUtilities::setLayoutSpacingAndMargins(layout, 0, 0);
    layout->addWidget(m_stackedWidget);
    layout->addStretch();
}

/**
 * Destructor.
 */
BrainBrowserWindowToolBarChartTwoAttributes::~BrainBrowserWindowToolBarChartTwoAttributes()
{
}

/**
 * Update content of this tool bar component.
 *
 * @param browserTabContent
 *     Content of the browser tab.
 */
void
BrainBrowserWindowToolBarChartTwoAttributes::updateContent(BrowserTabContent* browserTabContent)
{
    if (browserTabContent != NULL) {
        ModelChartTwo* modelChartTwo = browserTabContent->getDisplayedChartTwoModel();
        const int32_t tabIndex = browserTabContent->getTabNumber();
        if (modelChartTwo != NULL) {
            switch (modelChartTwo->getSelectedChartTwoDataType(tabIndex)) {
                case ChartTwoDataTypeEnum::CHART_DATA_TYPE_INVALID:
                    break;
                case ChartTwoDataTypeEnum::CHART_DATA_TYPE_HISTOGRAM:
                    break;
                case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_LAYER:
                    m_stackedWidget->setCurrentWidget(m_cartesianChartAttributesWidget);
                    m_cartesianChartAttributesWidget->updateContent();
                    break;
                case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES:
                    m_stackedWidget->setCurrentWidget(m_cartesianChartAttributesWidget);
                    m_cartesianChartAttributesWidget->updateContent();
                    break;
                case ChartTwoDataTypeEnum::CHART_DATA_TYPE_MATRIX:
                    m_stackedWidget->setCurrentWidget(m_matrixChartTwoAttributesWidget);
                    m_matrixChartTwoAttributesWidget->updateContent();
                    break;
            }
        }
    }
}

/**
 * @return Matrix chart interface in this widget.  Returned value will
 * be NULL if the chart (or no chart) is not a Matrix Chart.
 */
ChartTwoMatrixDisplayProperties*
BrainBrowserWindowToolBarChartTwoAttributes::getChartableTwoMatrixDisplayProperties()
{
    ChartTwoMatrixDisplayProperties* matrixDisplayProperties = NULL;

    BrowserTabContent* browserTabContent = getTabContentFromSelectedTab();
    if (browserTabContent != NULL) {
        matrixDisplayProperties = browserTabContent->getChartTwoMatrixDisplayProperties();
    }
    
    return matrixDisplayProperties;
}


/**
 * Update the graphics.
 */
void
BrainBrowserWindowToolBarChartTwoAttributes::updateGraphics()
{
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
}




/* ===========================================================================*/

/**
 * \class caret::CartesianChartTwoAttributesWidget
 * \brief Controls for cartesian chart attributes.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param brainBrowserWindowToolBarChartAttributes
 *   The parent attributes widget.
 */
CartesianChartTwoAttributesWidget::CartesianChartTwoAttributesWidget(BrainBrowserWindowToolBarChartTwoAttributes* brainBrowserWindowToolBarChartAttributes,
                                                                     const QString& parentObjectName)
: QWidget(brainBrowserWindowToolBarChartAttributes)
{
    m_brainBrowserWindowToolBarChartAttributes = brainBrowserWindowToolBarChartAttributes;
    
    QLabel* cartesianLineWidthLabel = new QLabel("Line width ");
    m_cartesianLineWidthDoubleSpinBox = WuQFactory::newDoubleSpinBoxWithMinMaxStepDecimalsSignalDouble(0.0,
                                                                                                       10000.0,
                                                                                                       0.1,
                                                                                                       1,
                                                                                                       this,
                                                                                                       SLOT(cartesianLineWidthValueChanged(double)));
    m_cartesianLineWidthDoubleSpinBox->setFixedWidth(65);
    m_cartesianLineWidthDoubleSpinBox->setToolTip("Width of line");
    m_cartesianLineWidthDoubleSpinBox->setObjectName(parentObjectName
                                                     + ":LineWidth");
    WuQMacroManager::instance()->addMacroSupportToObject(m_cartesianLineWidthDoubleSpinBox,
                                                         "Set cartesian chart line width");
    
    
    QGridLayout* gridLayout = new QGridLayout(this);
    WuQtUtilities::setLayoutSpacingAndMargins(gridLayout, 0, 0);
    gridLayout->setColumnStretch(0, 0);
    gridLayout->setColumnStretch(0, 100);
    gridLayout->addWidget(cartesianLineWidthLabel, 0, 0);
    gridLayout->addWidget(m_cartesianLineWidthDoubleSpinBox, 0, 1);
    
    this->setFixedSize(this->sizeHint());
}

/**
 * Destructor.
 */
CartesianChartTwoAttributesWidget::~CartesianChartTwoAttributesWidget()
{
    
}

/**
 * Update the content of this widget.
 */
void
CartesianChartTwoAttributesWidget::updateContent()
{
}

/**
 * Called when the cartesian line width is changed.
 *
 * @param value
 *    New value for line width.
 */
void
CartesianChartTwoAttributesWidget::cartesianLineWidthValueChanged(double /*value*/)
{
}




/* ===========================================================================*/



/**
 * \class caret::MatrixChartTwoAttributesWidget
 * \brief Controls for matrix chart attributes.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param brainBrowserWindowToolBarChartAttributes
 *   The parent attributes widget.
 */
MatrixChartTwoAttributesWidget::MatrixChartTwoAttributesWidget(BrainBrowserWindowToolBarChartTwoAttributes* brainBrowserWindowToolBarChartAttributes,
                                                               const QString& parentObjectName)
: QWidget(brainBrowserWindowToolBarChartAttributes),
EventListenerInterface()
{
    m_brainBrowserWindowToolBarChartAttributes = brainBrowserWindowToolBarChartAttributes;
        
    m_highlightSelectionCheckBox = new QCheckBox("Highlight Selection");
    m_highlightSelectionCheckBox->setToolTip("Highlight selected row/column in the matrix");
    QObject::connect(m_highlightSelectionCheckBox, SIGNAL(clicked(bool)),
                     this, SLOT(valueChanged()));
    m_highlightSelectionCheckBox->setObjectName(parentObjectName
                                                     + ":HighlightSelection");
    WuQMacroManager::instance()->addMacroSupportToObject(m_highlightSelectionCheckBox,
                                                         "Enable outline of selected chart matrix row/column");
    
    m_displayGridLinesCheckBox = new QCheckBox("Show Grid Outline");
    QObject::connect(m_displayGridLinesCheckBox, SIGNAL(clicked(bool)),
                     this, SLOT(valueChanged()));
    m_displayGridLinesCheckBox->setToolTip("Outline cells in the matrix");
    m_displayGridLinesCheckBox->setObjectName(parentObjectName
                                                     + ":EnableGridOutline");
    WuQMacroManager::instance()->addMacroSupportToObject(m_displayGridLinesCheckBox,
                                                         "Enable matrix chart grid outline");
    
    const int32_t COLUMN_LABEL  = 0;
    
    QGridLayout* gridLayout = new QGridLayout(this);
    WuQtUtilities::setLayoutSpacingAndMargins(gridLayout, 2, 2);
    int32_t rowIndex = gridLayout->rowCount();
    gridLayout->addWidget(m_highlightSelectionCheckBox, rowIndex, COLUMN_LABEL, 1, 2, Qt::AlignLeft);
    rowIndex++;
    gridLayout->addWidget(m_displayGridLinesCheckBox, rowIndex, COLUMN_LABEL, 1, 2, Qt::AlignLeft);
    rowIndex++;
    
    this->setFixedSize(this->sizeHint());
}

/**
 * Destructor.
 */
MatrixChartTwoAttributesWidget::~MatrixChartTwoAttributesWidget()
{
    EventManager::get()->removeAllEventsFromListener(this);
}

/**
 * Receive an event.
 *
 * @param event
 *    The event.
 */
void
MatrixChartTwoAttributesWidget::receiveEvent(Event* /*event*/)
{
}


/**
 * Update the content of this widget.
 */
void
MatrixChartTwoAttributesWidget::updateContent()
{
    ChartTwoMatrixDisplayProperties* matrixDisplayProperties = m_brainBrowserWindowToolBarChartAttributes->getChartableTwoMatrixDisplayProperties();
    if (matrixDisplayProperties != NULL) {
        m_highlightSelectionCheckBox->blockSignals(true);
        m_highlightSelectionCheckBox->setChecked(matrixDisplayProperties->isSelectedRowColumnHighlighted());
        m_highlightSelectionCheckBox->blockSignals(false);
        
        m_displayGridLinesCheckBox->blockSignals(true);
        m_displayGridLinesCheckBox->setChecked(matrixDisplayProperties->isGridLinesDisplayed());
        m_displayGridLinesCheckBox->blockSignals(false);
    }
}

/**
 * Gets called when user changes value of a user-interface component.
 */
void
MatrixChartTwoAttributesWidget::valueChanged()
{
    ChartTwoMatrixDisplayProperties* matrixDisplayProperties = m_brainBrowserWindowToolBarChartAttributes->getChartableTwoMatrixDisplayProperties();
    if (matrixDisplayProperties != NULL) {
        matrixDisplayProperties->setSelectedRowColumnHighlighted(m_highlightSelectionCheckBox->isChecked());
        matrixDisplayProperties->setGridLinesDisplayed(m_displayGridLinesCheckBox->isChecked());

        BrowserTabContent* tabContent = m_brainBrowserWindowToolBarChartAttributes->getTabContentFromSelectedTab();
        CaretAssert(tabContent);
        
        tabContent->updateChartModelYokedBrowserTabs();
        m_brainBrowserWindowToolBarChartAttributes->updateGraphics();
        EventManager::get()->sendEvent(EventUserInterfaceUpdate().addToolBar().getPointer());
    }
}


