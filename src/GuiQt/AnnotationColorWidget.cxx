
/*LICENSE_START*/
/*
 *  Copyright (C) 2015 Washington University School of Medicine
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

#define __ANNOTATION_COLOR_WIDGET_DECLARE__
#include "AnnotationColorWidget.h"
#undef __ANNOTATION_COLOR_WIDGET_DECLARE__

#include <QAction>
#include <QColorDialog>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QLabel>
#include <QToolButton>
#include <QVBoxLayout>

#include "AnnotationManager.h"
#include "AnnotationTwoCoordinateShape.h"
#include "AnnotationOneCoordinateShape.h"
#include "AnnotationRedoUndoCommand.h"
#include "Brain.h"
#include "BrainOpenGL.h"
#include "CaretAssert.h"
#include "CaretColorEnumMenu.h"
#include "EventAnnotationGetBeingDrawnInWindow.h"
#include "EventGraphicsPaintSoonAllWindows.h"
#include "EventManager.h"
#include "GuiManager.h"
#include "WuQFactory.h"
#include "WuQMessageBox.h"
#include "WuQWidgetObjectGroup.h"
#include "WuQtUtilities.h"

using namespace caret;


    
/**
 * \class caret::AnnotationColorWidget 
 * \brief Widget for annotation color selection.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param userInputMode
 *    The input mode
 * @param parentWidgetType
 *   Type of parent widget
 * @param browserWindowIndex
 *    Index of browser window
 * @param parent

 *     Parent for this widget.
 */
AnnotationColorWidget::AnnotationColorWidget(const UserInputModeEnum::Enum userInputMode,
                                             const AnnotationWidgetParentEnum::Enum parentWidgetType,
                                             const int32_t browserWindowIndex,
                                             QWidget* parent)
: QWidget(parent),
m_userInputMode(userInputMode),
m_parentWidgetType(parentWidgetType),
m_browserWindowIndex(browserWindowIndex)
{
    bool showFillFlag = true;
    switch (userInputMode) {
        case UserInputModeEnum::Enum::ANNOTATIONS:
            showFillFlag = true;
            break;
        case UserInputModeEnum::Enum::BORDERS:
            break;
        case UserInputModeEnum::Enum::FOCI:
            break;
        case UserInputModeEnum::Enum::IMAGE:
            break;
        case UserInputModeEnum::Enum::INVALID:
            break;
        case UserInputModeEnum::Enum::SAMPLES_EDITING:
            showFillFlag = false;
            break;
        case UserInputModeEnum::Enum::TILE_TABS_LAYOUT_EDITING:
            break;
        case UserInputModeEnum::Enum::VIEW:
            break;
        case UserInputModeEnum::Enum::VOLUME_EDIT:
            break;
    }
    
    QLabel* backFillLabel(new QLabel("Fill\nColor"));
    backFillLabel->setAlignment(Qt::AlignHCenter);
    
    QLabel* foreLineLabel      = new QLabel("Line\nColor");
    foreLineLabel->setAlignment(Qt::AlignHCenter);
    
    QLabel* lineLabel = new QLabel("Line\nWidth");
    lineLabel->setAlignment(Qt::AlignHCenter);
    
    const QSize toolButtonSize(16, 16);
    
    /*
     * Background color menu
     */
    m_backgroundColorMenu = new CaretColorEnumMenu((CaretColorEnum::OPTION_INCLUDE_CUSTOM_COLOR
                                                    | CaretColorEnum::OPTION_INCLUDE_NONE_COLOR));
    QObject::connect(m_backgroundColorMenu, SIGNAL(colorSelected(const CaretColorEnum::Enum)),
                     this, SLOT(backgroundColorSelected(const CaretColorEnum::Enum)));
    
    /*
     * Background action and tool button
     */
    m_backgroundColorAction = new QAction("B",
                                          this);
    m_backgroundColorAction->setToolTip("Adjust the fill color");
    m_backgroundColorAction->setMenu(m_backgroundColorMenu);
    m_backgroundToolButton = new QToolButton();
    m_backgroundToolButton->setDefaultAction(m_backgroundColorAction);
    m_backgroundToolButton->setIconSize(toolButtonSize);
    WuQtUtilities::setToolButtonStyleForQt5Mac(m_backgroundToolButton);
    
    /*
     * Widget/object group for background widgets
     */
    m_backgroundColorWidgetGroup = new WuQWidgetObjectGroup(this);
    m_backgroundColorWidgetGroup->add(backFillLabel);
    m_backgroundColorWidgetGroup->add(m_backgroundToolButton);

    /*
     * Line color menu
     */
    m_lineColorMenu = new CaretColorEnumMenu((CaretColorEnum::OPTION_INCLUDE_CUSTOM_COLOR
                                                    | CaretColorEnum::OPTION_INCLUDE_NONE_COLOR));
    QObject::connect(m_lineColorMenu, SIGNAL(colorSelected(const CaretColorEnum::Enum)),
                     this, SLOT(lineColorSelected(const CaretColorEnum::Enum)));
    
    /*
     * Line color action and toolbutton
     */
    m_lineColorAction = new QAction("F",
                                          this);
    m_lineColorAction->setToolTip("Adjust the line color");
    m_lineColorAction->setMenu(m_lineColorMenu);
    m_lineToolButton = new QToolButton();
    m_lineToolButton->setDefaultAction(m_lineColorAction);
    m_lineToolButton->setIconSize(toolButtonSize);
    WuQtUtilities::setToolButtonStyleForQt5Mac(m_lineToolButton);
    
    m_lineColorWidgetGroup = new WuQWidgetObjectGroup(this);
    m_lineColorWidgetGroup->add(foreLineLabel);
    m_lineColorWidgetGroup->add(m_lineColorMenu);
    m_lineColorWidgetGroup->add(m_lineToolButton);
    
    /*
     * Line thickness
     */
    float minimumLineWidthPercentage = 0.1;
    float maximumLineWidthPercentage = 100.0;
    float lineWidthStep = 0.1;
    m_lineThicknessWidgetGroup = new WuQWidgetObjectGroup(this);
    m_lineThicknessSpinBox = WuQFactory::newDoubleSpinBoxWithMinMaxStepDecimalsSignalDouble(minimumLineWidthPercentage,
                                                                                            maximumLineWidthPercentage,
                                                                                            lineWidthStep,
                                                                                            1,
                                                                                            this,
                                                                                            SLOT(lineThicknessSpinBoxValueChanged(double)));
    WuQtUtilities::setWordWrappedToolTip(m_lineThicknessSpinBox,
                                         "Adjust the line thickness as percentage of tab/window height");
    m_lineThicknessSpinBox->setSuffix("%");
    
    m_lineThicknessWidgetGroup->add(lineLabel);
    m_lineThicknessWidgetGroup->add(m_lineThicknessSpinBox);

    /*
     * Layout
     */
    QGridLayout* gridLayout = new QGridLayout(this);
    WuQtUtilities::setLayoutSpacingAndMargins(gridLayout, 2, 0);
    gridLayout->setVerticalSpacing(0);
    gridLayout->setHorizontalSpacing(2);
    
    switch (m_parentWidgetType) {
        case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
        {
            gridLayout->addWidget(lineLabel,
                                  0, 0,
                                  Qt::AlignHCenter);
            gridLayout->addWidget(foreLineLabel,
                                  0, 1,
                                  Qt::AlignHCenter);
            gridLayout->addWidget(m_lineThicknessSpinBox,
                                  2, 0,
                                  Qt::AlignHCenter);
            gridLayout->addWidget(m_lineToolButton,
                                  2, 1,
                                  Qt::AlignHCenter);
            if (showFillFlag) {
                gridLayout->addWidget(backFillLabel,
                                      0, 2,
                                      Qt::AlignHCenter);
                gridLayout->addWidget(m_backgroundToolButton,
                                      2, 2,
                                      Qt::AlignHCenter);
            }
            else {
                backFillLabel->setHidden(true);
                m_backgroundToolButton->setHidden(true);
            }
        }
            break;
        case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
            CaretAssert(0);
            break;
    }
    
    /*
     * Layout widgets
     */
    backgroundColorSelected(CaretColorEnum::WHITE);
    lineColorSelected(CaretColorEnum::BLACK);
    
    setSizePolicy(QSizePolicy::Fixed,
                  QSizePolicy::Fixed);
}

/**
 * Destructor.
 */
AnnotationColorWidget::~AnnotationColorWidget()
{
}

/**
 * Receive an event.
 *
 * @param event
 *     The event that the receive can respond to.
 */
void
AnnotationColorWidget::receiveEvent(Event* /*event*/)
{
}

/**
 * Update with the given annotation.
 *
 * @param annotationsIn
 */
void
AnnotationColorWidget::updateContent(std::vector<Annotation*>& annotationsIn)
{
    EventAnnotationGetBeingDrawnInWindow annDrawEvent(m_userInputMode,
                                                      m_browserWindowIndex);
    EventManager::get()->sendEvent(annDrawEvent.getPointer());
    m_annotationBeingDrawn = annDrawEvent.getAnnotation();

    std::vector<Annotation*> annotations(annotationsIn);
    if (annotations.empty()) {
        if (m_annotationBeingDrawn != NULL) {
            annotations.push_back(m_annotationBeingDrawn);
        }
    }
    else {
        m_annotationBeingDrawn = NULL;
    }

    m_lineColorAnnotations.clear();
    m_lineColorAnnotations.reserve(annotations.size());
    
    m_lineThicknessAnnotations.clear();
    m_lineThicknessAnnotations.reserve(annotations.size());

    m_backgroundColorAnnotations.clear();
    m_backgroundColorAnnotations.reserve(annotations.size());
    
    bool haveAnnotationsFlag = false;
    for (auto a : annotations) {
        if (a->testProperty(Annotation::Property::FILL_COLOR)) {
            m_backgroundColorAnnotations.push_back(a);
        }
        if (a->testProperty(Annotation::Property::LINE_COLOR)) {
            m_lineColorAnnotations.push_back(a);
        }
        if (a->testProperty(Annotation::Property::LINE_THICKNESS)) {
            m_lineThicknessAnnotations.push_back(a);
        }
        
        haveAnnotationsFlag = true;
    }
    
    setEnabled(haveAnnotationsFlag);
    
    updateBackgroundColorButton();
    updateLineColorButton();
    updateLineThicknessSpinBox();
}

/**
 * Gets called when the background color is changed.
 *
 * @param caretColor
 *     Color that was selected.
 */
void
AnnotationColorWidget::backgroundColorSelected(const CaretColorEnum::Enum caretColor)
{
    if (! m_backgroundColorAnnotations.empty()) {
        float rgba[4];
        m_backgroundColorAnnotations[0]->getCustomBackgroundColor(rgba);
        
        if (caretColor == CaretColorEnum::CUSTOM) {
            QColor color;
            color.setRgbF(rgba[0], rgba[1], rgba[2]);
            
            QColor newColor = QColorDialog::getColor(color,
                                                     m_backgroundToolButton,
                                                     "Background Color");
            if (newColor.isValid()) {
                rgba[0] = newColor.redF();
                rgba[1] = newColor.greenF();
                rgba[2] = newColor.blueF();
                

                switch (m_parentWidgetType) {
                    case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
                        Annotation::setUserDefaultCustomBackgroundColor(rgba);
                        break;
                    case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
                        CaretAssert(0);
                        break;
                }
            }
        }
        
        if ( ! isBothColorsSetToNoneAllowed(this,
                                            caretColor,
                                            m_lineColorMenu->getSelectedColor(),
                                            m_backgroundColorAnnotations)) {
            return;
        }
        
        if (m_annotationBeingDrawn) {
            m_annotationBeingDrawn->setBackgroundColor(caretColor);
            m_annotationBeingDrawn->setCustomBackgroundColor(rgba);
        }
        else {
            AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
            undoCommand->setModeColorBackground(caretColor,
                                                rgba,
                                                m_backgroundColorAnnotations);
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager(m_userInputMode);
            
            AString errorMessage;
            if ( ! annMan->applyCommand(undoCommand,
                                        errorMessage)) {
                WuQMessageBox::errorOk(this,
                                       errorMessage);
            }
        }

        Annotation::setUserDefaultBackgroundColor(caretColor);
        EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    }

    updateBackgroundColorButton();
    EventManager::get()->sendEvent(EventGraphicsPaintSoonAllWindows().getPointer());
}

/**
 * Update the background color.
 */
void
AnnotationColorWidget::updateBackgroundColorButton()
{
    CaretColorEnum::Enum colorEnum = CaretColorEnum::NONE;
    float rgba[4];
    CaretColorEnum::toRGBAFloat(colorEnum, rgba);
    rgba[3] = 1.0;
    
    bool enableBackgroundFlag = false;
    const int32_t numAnnotations = static_cast<int32_t>(m_backgroundColorAnnotations.size());
    if (numAnnotations > 0) {
        bool firstColorSupportFlag = true;
        bool allSameColorFlag = true;
        
        for (int32_t i = 0; i < numAnnotations; i++) {
            if (firstColorSupportFlag) {
                m_backgroundColorAnnotations[i]->getBackgroundColorRGBA(rgba);
                firstColorSupportFlag = false;
                enableBackgroundFlag = true;
            }
            else {
                float colorRGBA[4];
                m_backgroundColorAnnotations[i]->getBackgroundColorRGBA(colorRGBA);
                for (int32_t iColor = 0; iColor < 4; iColor++) {
                    if (rgba[iColor] != colorRGBA[iColor]) {
                        allSameColorFlag = false;
                        break;
                    }
                }
                
                if ( ! allSameColorFlag) {
                    break;
                }
            }
        }
        
        if (allSameColorFlag) {
            colorEnum = m_backgroundColorAnnotations[0]->getBackgroundColor();
            m_backgroundColorAnnotations[0]->getBackgroundColorRGBA(rgba);
            
            float customRGBA[4];
            m_backgroundColorAnnotations[0]->getCustomBackgroundColor(customRGBA);
            m_backgroundColorMenu->setCustomIconColor(customRGBA);
            
            switch (m_parentWidgetType) {
                case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
                    Annotation::setUserDefaultBackgroundColor(colorEnum);
                    Annotation::setUserDefaultCustomBackgroundColor(customRGBA);
                    break;
                case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
                    CaretAssert(0);
                    break;
            }
        }
        
        
    }
    
    m_backgroundColorWidgetGroup->setEnabled(enableBackgroundFlag);
    if ( ! enableBackgroundFlag) {
        colorEnum = CaretColorEnum::NONE;
    }
    
    QPixmap pm = WuQtUtilities::createCaretColorEnumPixmap(m_backgroundToolButton,
                                                           24, 24,
                                                           colorEnum,
                                                           rgba,
                                                           false);
    QIcon icon(pm);
    
    m_backgroundColorAction->setIcon(icon);
    m_backgroundColorMenu->setSelectedColor(colorEnum);
    
    /*
     * Used for event requesting color for drawing a new annotation
     */
    m_currentBackgroundColorForDrawingNewAnnotation.setCaretColorEnum(colorEnum);
    std::array<uint8_t, 4> rgbaByte {
        static_cast<uint8_t>(rgba[0] * 255.0),
        static_cast<uint8_t>(rgba[1] * 255.0),
        static_cast<uint8_t>(rgba[2] * 255.0),
        255
    };
    m_currentBackgroundColorForDrawingNewAnnotation.setCustomColorRGBA(rgbaByte);
}


/**
 * Update the line color.
 */
void
AnnotationColorWidget::updateLineColorButton()
{
    CaretColorEnum::Enum colorEnum = CaretColorEnum::NONE;
    float rgba[4];
    CaretColorEnum::toRGBAFloat(colorEnum, rgba);
    rgba[3] = 1.0;
    
    const int32_t numAnnotations = static_cast<int32_t>(m_lineColorAnnotations.size());
    bool enableLineFlag = false;
    if (numAnnotations > 0) {
        bool firstColorSupportFlag = true;
        bool allSameColorFlag = true;
        
        for (int32_t i = 0; i < numAnnotations; i++) {
                if (firstColorSupportFlag) {
                    m_lineColorAnnotations[i]->getLineColorRGBA(rgba);
                    firstColorSupportFlag = false;
                    enableLineFlag = true;
                }
                else {
                    float colorRGBA[4];
                    m_lineColorAnnotations[i]->getLineColorRGBA(colorRGBA);
                    for (int32_t iColor = 0; iColor < 4; iColor++) {
                        if (rgba[iColor] != colorRGBA[iColor]) {
                            allSameColorFlag = false;
                            break;
                        }
                    }
                    
                    if ( ! allSameColorFlag) {
                        break;
                    }
                }
        }
        
        if (allSameColorFlag) {
            colorEnum = m_lineColorAnnotations[0]->getLineColor();
            m_lineColorAnnotations[0]->getLineColorRGBA(rgba);
            
            float customRGBA[4];
            m_lineColorAnnotations[0]->getCustomLineColor(customRGBA);
            m_lineColorMenu->setCustomIconColor(customRGBA);

            switch (m_parentWidgetType) {
                case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
                    setUserDefaultLineColor(colorEnum,
                                            customRGBA);
                    break;
                case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
                    CaretAssert(0);
                    break;
            }
            
        }
    }

    m_lineColorWidgetGroup->setEnabled(enableLineFlag);
    
    if ( ! enableLineFlag) {
        colorEnum = CaretColorEnum::NONE;
    }
    
    QPixmap pm = WuQtUtilities::createCaretColorEnumPixmap(m_lineToolButton, 24, 24, colorEnum, rgba, true);
    m_lineColorAction->setIcon(QIcon(pm));
    m_lineColorMenu->setSelectedColor(colorEnum);
    
    /*
     * Used for event requesting color for drawing a new annotation
     */
    m_currentLineColorForDrawingNewAnnotation.setCaretColorEnum(colorEnum);
    std::array<uint8_t, 4> rgbaByte {
        static_cast<uint8_t>(rgba[0] * 255.0),
        static_cast<uint8_t>(rgba[1] * 255.0),
        static_cast<uint8_t>(rgba[2] * 255.0),
        255
    };
    m_currentLineColorForDrawingNewAnnotation.setCustomColorRGBA(rgbaByte);
}

/**
 * Tests for both background and line colors both set to none.  This is allowed
 * for some annotations types such as text, which uses a text color, or images.
 *
 * @param widget
 *     Widget on which error dialog is displayed.
 * @param colorOne
 *     One of the background or line colors.
 * @param colorTwo
 *     The other one of the background or line colors.
 * @param annotations
 *     The selected annotations.
 * @return
 *     True if the colors are acceptable, else false (and a an error
 *     message dialog is displayed).
 */
bool
AnnotationColorWidget::isBothColorsSetToNoneAllowed(QWidget* widget,
                                                    const CaretColorEnum::Enum colorOne,
                                                    const CaretColorEnum::Enum colorTwo,
                                                    const std::vector<Annotation*>& annotations)
{
    if ((colorOne == CaretColorEnum::NONE)
        && (colorTwo == CaretColorEnum::NONE)) {
        
        bool allowBothColorsNoneFlag = true;
        
        for (std::vector<Annotation*>::const_iterator iter = annotations.begin();
             iter != annotations.end();
             iter++) {
            const Annotation* ann = *iter;
            CaretAssert(ann);
            switch (ann->getType()) {
                case AnnotationTypeEnum::BOX:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::BROWSER_TAB:
                    break;
                case AnnotationTypeEnum::COLOR_BAR:
                    break;
                case AnnotationTypeEnum::IMAGE:
                    break;
                case AnnotationTypeEnum::LINE:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::OVAL:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::POLYHEDRON:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::POLYGON:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::POLYLINE:
                    allowBothColorsNoneFlag = false;
                    break;
                case AnnotationTypeEnum::SCALE_BAR:
                    break;
                case AnnotationTypeEnum::TEXT:
                    break;
            }
        }
        
        if ( ! allowBothColorsNoneFlag) {
            const AString message("Setting both Line and Fill colors to NONE is not allowed for the selected annotation(s).");
            WuQMessageBox::errorOk(widget,
                                   message);
            return false;
        }
    }

    return true;
}


/**
 * Gets called when the line color is changed.
 *
 * @param caretColor
 *     Color that was selected.
 */
void
AnnotationColorWidget::lineColorSelected(const CaretColorEnum::Enum caretColor)
{
    if ( ! m_lineColorAnnotations.empty()) {
        CaretAssertVectorIndex(m_lineColorAnnotations, 0);
        
        float rgba[4];
        m_lineColorAnnotations[0]->getCustomLineColor(rgba);
        
        if (caretColor == CaretColorEnum::CUSTOM) {
            QColor color;
            color.setRgbF(rgba[0], rgba[1], rgba[2]);
            
            QColor newColor = QColorDialog::getColor(color,
                                                     m_backgroundToolButton,
                                                     "Line Color");
            if (newColor.isValid()) {
                rgba[0] = newColor.redF();
                rgba[1] = newColor.greenF();
                rgba[2] = newColor.blueF();
            }
        }
        
        if ( ! isBothColorsSetToNoneAllowed(this,
                                            caretColor,
                                            m_backgroundColorMenu->getSelectedColor(),
                                            m_lineColorAnnotations)) {
            return;
        }
        
        if (m_annotationBeingDrawn != NULL) {
            m_annotationBeingDrawn->setLineColor(caretColor);
            m_annotationBeingDrawn->setCustomLineColor(rgba);
        }
        else {
            AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
            undoCommand->setModeColorLine(caretColor,
                                          rgba,
                                          m_lineColorAnnotations);
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager(m_userInputMode);
            
            AString errorMessage;
            if ( ! annMan->applyCommand(undoCommand,
                                        errorMessage)) {
                WuQMessageBox::errorOk(this,
                                       errorMessage);
            }
        }
        
        setUserDefaultLineColor(caretColor,
                                rgba);
    }
    
    updateLineColorButton();
    
    switch (m_parentWidgetType) {
        case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
            EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
            break;
        case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
            break;
    }
    
    EventManager::get()->sendEvent(EventGraphicsPaintSoonAllWindows().getPointer());
}

/**
 * Set the user default line color.
 *
 * @param caretColor
 *     The line color.
 * @param customRGBA
 *     The custom line color RGBA components.
 */
void
AnnotationColorWidget::setUserDefaultLineColor(const CaretColorEnum::Enum caretColor,
                                               const float customRGBA[4])
{
    if ( ! m_lineColorAnnotations.empty()) {
        bool allTextFlag  = true;
        bool someTextFlag = false;
        
        for (std::vector<Annotation*>::iterator iter = m_lineColorAnnotations.begin();
             iter != m_lineColorAnnotations.end();
             iter++) {
            const Annotation* ann = *iter;
            if (ann->getType() == AnnotationTypeEnum::TEXT) {
                someTextFlag = true;
            }
            else {
                allTextFlag = false;
            }
        }
        
        /*
         * Note: Text has its own default line color.  Without it, if the
         * user creates a text annotation it will get a box around it since
         * other annotations frequently use a line color.
         */
        if (allTextFlag
            || someTextFlag) {
            Annotation::setUserDefaultForTextLineColor(caretColor);
            if (caretColor == CaretColorEnum::CUSTOM) {
                Annotation::setUserDefaultForTextCustomLineColor(customRGBA);
            }
        }
        
        if (! allTextFlag) {
            Annotation::setUserDefaultLineColor(caretColor);
            if (caretColor == CaretColorEnum::CUSTOM) {
                Annotation::setUserDefaultCustomLineColor(customRGBA);
            }
        }
    }
}


/**
 * Gets called when the line thickness value changes.
 *
 * @param value
 *     New value for line thickness.
 */
void
AnnotationColorWidget::lineThicknessSpinBoxValueChanged(double value)
{
    if ( ! m_lineThicknessSpinBox->specialValueText().isEmpty()) {
        if (m_lineThicknessSpinBox->specialValueText()
            == m_lineThicknessSpinBox->text()) {
            /*
             * Ignore special text which is available when 
             * there are multiple annotations with different
             * line thicknesses.
             */
            std::cout << "Ignoring special text " << std::endl;
            return;
        }
    }
    
    if (m_annotationBeingDrawn != NULL) {
        m_annotationBeingDrawn->setLineWidthPercentage(value);
    }
    else {
        AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
        undoCommand->setModeLineWidth(value,
                                      m_lineThicknessAnnotations);
        AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager(m_userInputMode);
        
        AString errorMessage;
        if ( ! annMan->applyCommand(undoCommand,
                                    errorMessage)) {
            WuQMessageBox::errorOk(this,
                                   errorMessage);
        }
    }

    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    Annotation::setUserDefaultLineWidthPercentage(value);
    
    EventManager::get()->sendEvent(EventGraphicsPaintSoonAllWindows().getPointer());
}

/**
 * Update the line thickness spin box.
 */
void
AnnotationColorWidget::updateLineThicknessSpinBox()
{
    if (m_lineThicknessSpinBox == NULL) {
        return;
    }
    
    float lineWidthValue = 1.0;
    bool  lineWidthValid = false;
    bool  haveMultipleLineWidthValues = false;
    
    const int32_t numAnnotations = static_cast<int32_t>(m_lineThicknessAnnotations.size());
    if (numAnnotations > 0) {
        for (int32_t i = 0; i < numAnnotations; i++) {
            const float annLineWidth = m_lineThicknessAnnotations[i]->getLineWidthPercentage();
            if (lineWidthValid) {
                if (annLineWidth != lineWidthValue) {
                    haveMultipleLineWidthValues = true;
                }
                lineWidthValue = std::min(lineWidthValue,
                                          annLineWidth);
            }
            else {
                lineWidthValue = annLineWidth;
                lineWidthValid = true;
            }
        }
        
        if (lineWidthValid) {
            switch (m_parentWidgetType) {
                case AnnotationWidgetParentEnum::ANNOTATION_TOOL_BAR_WIDGET:
                    Annotation::setUserDefaultLineWidthPercentage(lineWidthValue);
                    break;
                case AnnotationWidgetParentEnum::PARENT_ENUM_FOR_LATER_USE:
                    CaretAssert(0);
                    break;
            }
        }
    }
    
    /*
     * When the selected annotations have different line
     * widths, the valid displayed is the minimum line
     * width with a suffix consisting of a plus symbol.
     */
    m_lineThicknessSpinBox->blockSignals(true);
    m_lineThicknessSpinBox->setValue(lineWidthValue);
    m_lineThicknessSpinBox->setEnabled(lineWidthValid);
    if (haveMultipleLineWidthValues) {
        m_lineThicknessSpinBox->setSuffix("%+");
    }
    else {
        m_lineThicknessSpinBox->setSuffix("%");
    }
    m_lineThicknessSpinBox->blockSignals(false);
    
    m_lineThicknessWidgetGroup->setEnabled(numAnnotations > 0);
}

