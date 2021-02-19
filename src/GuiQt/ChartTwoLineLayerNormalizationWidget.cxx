
/*LICENSE_START*/
/*
 *  Copyright (C) 2020 Washington University School of Medicine
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

#define __CHART_TWO_LINE_LAYER_NORMALIZATION_WIDGET_DECLARE__
#include "ChartTwoLineLayerNormalizationWidget.h"
#undef __CHART_TWO_LINE_LAYER_NORMALIZATION_WIDGET_DECLARE__

#include <limits>

#include <QCheckBox>
#include <QGridLayout>
#include <QLabel>
#include <QDoubleSpinBox>

#include "CaretAssert.h"
#include "ChartTwoDataCartesian.h"
#include "ChartTwoOverlay.h"
#include "EventGraphicsUpdateAllWindows.h"
#include "EventUserInterfaceUpdate.h"
#include "EventManager.h"
#include "GraphicsPrimitiveV3f.h"
#include "WuQtUtilities.h"

using namespace caret;

/**
 * \class caret::ChartTwoLineLayerNormalizationWidget
 * \brief Widget for editing chart line layer normalization parameters
 * \ingroup GuiQt
 */

/**
 * Constructor.
 */
ChartTwoLineLayerNormalizationWidget::ChartTwoLineLayerNormalizationWidget()
: QWidget()
{
    m_newMeanEnabledCheckBox = new QCheckBox("New Mean");
    QObject::connect(m_newMeanEnabledCheckBox, &QCheckBox::clicked,
                     this, &ChartTwoLineLayerNormalizationWidget::newMeanEnabledCheckBoxClicked);
    m_newMeanSpinBox = new QDoubleSpinBox();
    m_newMeanSpinBox->setDecimals(4);
    m_newMeanSpinBox->setSingleStep(0.1);
    m_newMeanSpinBox->setRange(-std::numeric_limits<float>::max(),
                               std::numeric_limits<float>::max());
    m_newMeanSpinBox->setMaximumWidth(100);
    QObject::connect(m_newMeanSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                     this, &ChartTwoLineLayerNormalizationWidget::newMeanValueChanged);

    m_newDeviationEnabledCheckBox = new QCheckBox("New Deviation");
    QObject::connect(m_newDeviationEnabledCheckBox, &QCheckBox::clicked,
                     this, &ChartTwoLineLayerNormalizationWidget::newDeviationEnabledCheckBoxClicked);
    m_newDeviationSpinBox = new QDoubleSpinBox();
    m_newDeviationSpinBox->setDecimals(4);
    m_newDeviationSpinBox->setSingleStep(0.1);
    m_newDeviationSpinBox->setRange(-std::numeric_limits<float>::max(),
                                    std::numeric_limits<float>::max());
    m_newDeviationSpinBox->setMaximumWidth(100);
    QObject::connect(m_newDeviationSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                     this, &ChartTwoLineLayerNormalizationWidget::newDeviationValueChanged);

    m_absoluteValueEnabledCheckBox = new QCheckBox("Absolute Values");
    QObject::connect(m_absoluteValueEnabledCheckBox, &QCheckBox::clicked,
                     this, &ChartTwoLineLayerNormalizationWidget::absoluteValueEnabledCheckBoxClicked);

    m_meanDevLabel = new QLabel("");
    
    QLabel* descripionLabel = new QLabel(GraphicsPrimitive::getNewMeanDeviationOperationDescriptionInHtml());
    descripionLabel->setWordWrap(true);
    descripionLabel->setMaximumWidth(400);
    
    QGridLayout* layout = new QGridLayout(this);
    layout->setColumnStretch(0, 0);
    layout->setColumnStretch(1, 100);
    int32_t row(0);
    layout->addWidget(m_absoluteValueEnabledCheckBox, row, 0, 1, 2, Qt::AlignLeft);
    row++;
    layout->addWidget(m_newMeanEnabledCheckBox, row, 0);
    layout->addWidget(m_newMeanSpinBox, row, 1, Qt::AlignLeft);
    row++;
    layout->addWidget(m_newDeviationEnabledCheckBox, row, 0);
    layout->addWidget(m_newDeviationSpinBox, row, 1, Qt::AlignLeft);
    row++;
    layout->addWidget(m_meanDevLabel, row, 0, 1, 2, Qt::AlignLeft);
    row++;
    layout->addWidget(WuQtUtilities::createHorizontalLineWidget(), row, 0, 1, 2);
    row++;
    layout->addWidget(descripionLabel, row, 0, 1, 2, Qt::AlignLeft);
    row++;

    m_blockUpdatesFlag = false;
}

/**
 * Destructor.
 */
ChartTwoLineLayerNormalizationWidget::~ChartTwoLineLayerNormalizationWidget()
{
}

/**
 * Update the content of the widget
 * @param chartTwoOverlay
 */
void
ChartTwoLineLayerNormalizationWidget::updateContent(ChartTwoOverlay* chartTwoOverlay)
{
    if (m_blockUpdatesFlag) {
        return;
    }
    
    m_chartTwoOverlay = chartTwoOverlay;
    
    QString meanDevText;
    
    bool validFlag(false);
    if (m_chartTwoOverlay != NULL) {
        m_absoluteValueEnabledCheckBox->setChecked(m_chartTwoOverlay->isLineChartNormalizationAbsoluteValueEnabled());
        
        m_newMeanEnabledCheckBox->setChecked(m_chartTwoOverlay->isLineChartNewMeanEnabled());
        QSignalBlocker meanBlocker(m_newMeanSpinBox);
        m_newMeanSpinBox->setValue(m_chartTwoOverlay->getLineChartNewMeanValue());
        
        m_newDeviationEnabledCheckBox->setChecked(m_chartTwoOverlay->isLineChartNewDeviationEnabled());
        QSignalBlocker devBlocker(m_newDeviationSpinBox);
        m_newDeviationSpinBox->setValue(m_chartTwoOverlay->getLineChartNewDeviationValue());
        
        float mean(0.0), dev(0.0);
        const ChartTwoDataCartesian* cartData = m_chartTwoOverlay->getLineLayerChartMapFileCartesianData();
        CaretAssert(cartData);
        const GraphicsPrimitiveV3f* primitive = cartData->getGraphicsPrimitive();
        primitive->getMeanAndStandardDeviationForY(mean, dev);
        
        const QString muCharacter("Data Mean"); //QChar(0x03BC));
        const QString sigmaCharacter("Data Deviation"); //QChar(0x03C3));
        meanDevText = (muCharacter + ": "
                       + QString::number(mean, 'f', 4)
                       + ", " + sigmaCharacter + ": "
                       + QString::number(dev, 'f', 4));
        validFlag = true;
    }
    
    this->setEnabled(validFlag);

    m_meanDevLabel->setText(meanDevText);
}

/**
 * Called when a value is changed
 */
void
ChartTwoLineLayerNormalizationWidget::valueChanged()
{
    if (m_chartTwoOverlay != NULL) {
        m_chartTwoOverlay->setLineChartNewMeanEnabled(m_newMeanEnabledCheckBox->isChecked());
        m_chartTwoOverlay->setLineChartNewMeanValue(m_newMeanSpinBox->value());
        m_chartTwoOverlay->setLineChartNewDeviationEnabled(m_newDeviationEnabledCheckBox->isChecked());
        m_chartTwoOverlay->setLineChartNewDeviationValue(m_newDeviationSpinBox->value());
        
        updateGraphics();
    }
}

/**
 * Called when new mean checkbox is clicked
 * @param clicked
 *    New clicked status
 */
void
ChartTwoLineLayerNormalizationWidget::newMeanEnabledCheckBoxClicked(bool clicked)
{
    m_chartTwoOverlay->setLineChartNewMeanEnabled(clicked);
    updateGraphics();
    updateToolBarChartAxes();
}

/**
 * Called when new deviation checkbox is clicked
 * @param clicked
 *    New clicked status
 */
void
ChartTwoLineLayerNormalizationWidget::newDeviationEnabledCheckBoxClicked(bool clicked)
{
    m_chartTwoOverlay->setLineChartNewDeviationEnabled(clicked);
    updateGraphics();
    updateToolBarChartAxes();
}

/**
 * Called when new deviation checkbox is clicked
 * @param clicked
 *    New clicked status
 */
void
ChartTwoLineLayerNormalizationWidget::absoluteValueEnabledCheckBoxClicked(bool clicked)
{
    m_chartTwoOverlay->setLineChartNormalizationAbsoluteValueEnabled(clicked);
    updateGraphics();
    updateToolBarChartAxes();
}

/**
 * Called when new mean value is changed
 * @param value
 *    New value
 */
void
ChartTwoLineLayerNormalizationWidget::newMeanValueChanged(double value)
{
    m_chartTwoOverlay->setLineChartNewMeanValue(value);
    updateGraphics();
    updateToolBarChartAxes();
}

/**
 * Called when new deviation value is changed
 * @param value
 *    New value
 */
void
ChartTwoLineLayerNormalizationWidget::newDeviationValueChanged(double value)
{
    m_chartTwoOverlay->setLineChartNewDeviationValue(value);
    updateGraphics();
    updateToolBarChartAxes();
}

/**
 * Called when new mean value is changed
 * @param value
 *    New value
 */
void
ChartTwoLineLayerNormalizationWidget::updateGraphics()
{
    
    m_blockUpdatesFlag = true;
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    m_blockUpdatesFlag = false;
}

/**
 * Called when new mean value is changed
 * @param value
 *    New value
 */
void
ChartTwoLineLayerNormalizationWidget::updateToolBarChartAxes()
{
    m_blockUpdatesFlag = true;
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_TOOLBAR_CHART_ORIENTED_AXES_UPDATE);
    m_blockUpdatesFlag = false;
    
}



