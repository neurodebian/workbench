
/*LICENSE_START*/
/*
 *  Copyright (C) 2016 Washington University School of Medicine
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

#define __ANNOTATION_PASTE_DIALOG_DECLARE__
#include "AnnotationPasteDialog.h"
#undef __ANNOTATION_PASTE_DIALOG_DECLARE__

#include <cmath>

#include <QGroupBox>
#include <QLabel>
#include <QVBoxLayout>

#include "Annotation.h"
#include "AnnotationCoordinateInformation.h"
#include "AnnotationCoordinateSelectionWidget.h"
#include "AnnotationFile.h"
#include "AnnotationManager.h"
#include "AnnotationMultiCoordinateShape.h"
#include "AnnotationTwoCoordinateShape.h"
#include "AnnotationPercentSizeText.h"
#include "AnnotationRedoUndoCommand.h"
#include "Brain.h"
#include "BrainBrowserWindow.h"
#include "BrainOpenGLViewportContent.h"
#include "BrainOpenGLWidget.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "EventGraphicsUpdateAllWindows.h"
#include "EventUserInterfaceUpdate.h"
#include "EventManager.h"
#include "GuiManager.h"
#include "MathFunctions.h"
#include "ModelSurfaceMontage.h"
#include "MouseEvent.h"
#include "WuQMessageBox.h"
#include "WuQtUtilities.h"

using namespace caret;


    
/**
 * \class caret::AnnotationPasteDialog 
 * \brief Dialog for pasting annotations.
 * \ingroup GuiQt
 */

/**
 * Paste the annotation on the clipboard using the mouse information.
 *
 * @param mouseEvent
 *     Information about where to paste the annotation.
 * @param windowIndex
 *     Window in which annotation is pasted.
 * @return
 *     Pointer to annotation that was pasted or NULL if the annotation
 *     on the clipboard was not pasted.
 */
Annotation*
AnnotationPasteDialog::pasteAnnotationOnClipboard(const MouseEvent& mouseEvent,
                                                  const int32_t windowIndex)
{
    Annotation* newPastedAnnotation = NULL;
    
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    if (annotationManager->isAnnotationOnClipboardValid()) {
        AnnotationFile* annotationFile = annotationManager->getAnnotationFileOnClipboard();
        Annotation* annotation = annotationManager->getAnnotationOnClipboard()->clone();
        
        BrainOpenGLViewportContent* viewportContent = mouseEvent.getViewportContent();
        AnnotationCoordinateInformation coordInfo;
        
        AnnotationCoordinateInformation::createCoordinateInformationFromXY(mouseEvent.getOpenGLWidget(),
                                                                        viewportContent,
                                                                        mouseEvent.getX(),
                                                                        mouseEvent.getY(),
                                                                        coordInfo);
        
        bool validCoordsFlag = false;
        
        AnnotationTwoCoordinateShape* oneDimShape = dynamic_cast<AnnotationTwoCoordinateShape*>(annotation);
        if (oneDimShape != NULL) {
            /*
             * Pasting line while preserving its orientation only
             * works for tab and window spaces.
             */
            validCoordsFlag = pasteOneDimensionalShape(oneDimShape,
                                                       coordInfo);
        }
        
        AnnotationMultiCoordinateShape* multiCoordShape = annotation->castToMultiCoordinateShape();
        if (multiCoordShape != NULL) {
            /*
             * Pasting poly line while only works for tab and window spaces
             */
            validCoordsFlag = pasteMultiCoordinateShape(multiCoordShape,
                                                        coordInfo);
        }
        
        std::vector<std::unique_ptr<AnnotationCoordinateInformation>> emptyMultiCoordInfo;
        if (! validCoordsFlag) {
            validCoordsFlag = AnnotationCoordinateInformation::setAnnotationCoordinatesForSpace(annotation,
                                                                                                annotation->getCoordinateSpace(),
                                                                                                &coordInfo,
                                                                                                NULL,
                                                                                                emptyMultiCoordInfo);
        }
        
        
        if (validCoordsFlag) {
            AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
            undoCommand->setModePasteAnnotation(annotationFile,
                                                annotation);
            switch (annotation->getType()) {
                case AnnotationTypeEnum::BOX:
                    break;
                case AnnotationTypeEnum::BROWSER_TAB:
                    CaretAssert(0);
                    break;
                case AnnotationTypeEnum::COLOR_BAR:
                    break;
                case AnnotationTypeEnum::IMAGE:
                    break;
                case AnnotationTypeEnum::LINE:
                    break;
                case AnnotationTypeEnum::OVAL:
                    break;
                case AnnotationTypeEnum::POLY_LINE:
                    break;
                case AnnotationTypeEnum::SCALE_BAR:
                    break;
                case AnnotationTypeEnum::TEXT:
                    break;
            }
            AString errorMessage;
            if ( ! annotationManager->applyCommand(UserInputModeEnum::Enum::ANNOTATIONS,
                                                   undoCommand,
                                        errorMessage)) {
                WuQMessageBox::errorOk(mouseEvent.getOpenGLWidget(),
                                       errorMessage);
            }
            newPastedAnnotation = annotation;
            
            annotationManager->selectAnnotationForEditing(windowIndex,
                                                AnnotationManager::SELECTION_MODE_SINGLE,
                                                false,
                                                annotation);
        }
        else {
            /*
             * Pasting annotation in its coordinate failed (user may have tried to paste
             * an annotation in surface space where there is no surface).
             */
            delete annotation;
            annotation = NULL;
            
            const QString message("The location for pasting the annotation is incompatible with the "
                                  "coordinate space "
                                  "used by the annotation on the clipboard.  Choose one of the coordinate "
                                  "spaces below to paste the annotation or press Cancel to cancel pasting "
                                  "of the annotation.");
            
            /*
             * The annotation dialog will create a new annotation for pasting
             */
            AnnotationPasteDialog pasteDialog(mouseEvent,
                                              annotationFile,
                                              annotationManager->getAnnotationOnClipboard(),
                                              message,
                                              mouseEvent.getOpenGLWidget());
            if (pasteDialog.exec() == AnnotationPasteDialog::Accepted) {
                newPastedAnnotation = pasteDialog.getAnnotationThatWasCreated();
            }
        }
    }
    
    return newPastedAnnotation;
}

/**
 * Paste the annotation on the clipboard using the mouse information
 * and allow the user to change the coordinate space.
 *
 * @param mouseEvent
 *     Information about where to paste the annotation.
 * @return
 *     Pointer to annotation that was pasted or NULL if the annotation
 *     on the clipboard was not pasted.
 */
Annotation*
AnnotationPasteDialog::pasteAnnotationOnClipboardChangeSpace(const MouseEvent& mouseEvent)
{
    Annotation* newPastedAnnotation = NULL;
    
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    if (annotationManager->isAnnotationOnClipboardValid()) {
        AnnotationFile* annotationFile = annotationManager->getAnnotationFileOnClipboard();
        
        AString message("Choose one of the coordinate "
                        "spaces below to paste the annotation or press Cancel to cancel pasting "
                        "of the annotation.");
        AnnotationPasteDialog pasteDialog(mouseEvent,
                                          annotationFile,
                                          annotationManager->getAnnotationOnClipboard(),
                                          message,
                                          mouseEvent.getOpenGLWidget());
        if (pasteDialog.exec() == AnnotationPasteDialog::Accepted) {
            newPastedAnnotation = pasteDialog.getAnnotationThatWasCreated();
        }
    }
    
    return newPastedAnnotation;
}

/**
 * Constructor.
 *
 * @param mouseEvent
 *     Information about where mouse was clicked.
 * @param annotationFile
 *     File that contains the annotation.
 * @param annotation
 *     Annotation that is copied and pasted.
 * @param informationMessage
 *     Message shown on dialog.
 * @param parent
 *     Parent widget of dialog.
 */
AnnotationPasteDialog::AnnotationPasteDialog(const MouseEvent& mouseEvent,
                                             AnnotationFile* annotationFile,
                                             const Annotation* annotation,
                                             const AString& informationMessage,
                                             QWidget* parent)
: WuQDialogModal("Paste Annotation",
                 parent),
m_mouseEvent(mouseEvent),
m_annotationFile(annotationFile),
m_annotation(annotation),
m_annotationThatWasCreated(NULL)
{
    CaretAssert(m_annotationFile);
    CaretAssert(m_annotation);
    
    /*
     * Get coordinates at the mouse location.
     */
    AnnotationCoordinateInformation::createCoordinateInformationFromXY(mouseEvent,
                                                                    m_coordInfo);
    m_coordinateSelectionWidget = new AnnotationCoordinateSelectionWidget(m_annotation->getType(),
                                                                          m_coordInfo,
                                                                          NULL);
    m_coordinateSelectionWidget->selectCoordinateSpace(m_annotation->getCoordinateSpace());
    
    QLabel* messageLabel = new QLabel(informationMessage);
    messageLabel->setWordWrap(true);
    
    QLabel* spaceLabel = new QLabel("Space of Annotation on Clipboard: "
                                    + AnnotationCoordinateSpaceEnum::toGuiName(m_annotation->getCoordinateSpace()));
    spaceLabel->setWordWrap(false);
    
    QGroupBox* coordGroupBox = new QGroupBox("Coordinate Space");
    QVBoxLayout* coordGroupLayout = new QVBoxLayout(coordGroupBox);
    coordGroupLayout->setMargin(0);
    coordGroupLayout->addWidget(m_coordinateSelectionWidget);
    
    QWidget* dialogWidget = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(dialogWidget);
    layout->addWidget(spaceLabel);
    layout->addSpacing(10);
    layout->addWidget(messageLabel);
    layout->addSpacing(10);
    layout->addWidget(coordGroupBox);
    
    setCentralWidget(dialogWidget,
                     SCROLL_AREA_NEVER);
}

/**
 * Destructor.
 */
AnnotationPasteDialog::~AnnotationPasteDialog()
{
}

/**
 * @return Get the annotation that was created.
 */
Annotation*
AnnotationPasteDialog::getAnnotationThatWasCreated()
{
    return m_annotationThatWasCreated;
}

/**
 * Paste a one-dimensional shape (line) keeping its start to end
 * coordinate orientation.  The start coordinate is pasted at
 * the coordinate in 'coordInfo'.
 *
 * @param oneDimShape
 *     One dimensional shape that will be pasted.
 * @param coordInfo
 *     Coordinate information that will be used for the shape's 'start' coordinate.
 * @return
 *     True if the shape's coordinate was updated for pasting, else false.
 */
bool
AnnotationPasteDialog::pasteOneDimensionalShape(AnnotationTwoCoordinateShape* oneDimShape,
                                                   AnnotationCoordinateInformation& coordInfo)
{
    
    bool tabFlag = false;
    bool windowFlag = false;
    
    switch (oneDimShape->getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            tabFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            windowFlag = true;
            break;
    }
    
    bool validCoordsFlag = false;
    
    if (tabFlag
        || windowFlag) {
        float startXYZ[3];
        float endXYZ[3];
        oneDimShape->getStartCoordinate()->getXYZ(startXYZ);
        oneDimShape->getEndCoordinate()->getXYZ(endXYZ);
        const float diffXYZ[3] = {
            endXYZ[0] - startXYZ[0],
            endXYZ[1] - startXYZ[1],
            endXYZ[2] - startXYZ[2]
        };
        
        if (tabFlag
            && (coordInfo.m_tabSpaceInfo.m_index >= 0)) {
            startXYZ[0] = coordInfo.m_tabSpaceInfo.m_xyz[0];
            startXYZ[1] = coordInfo.m_tabSpaceInfo.m_xyz[1];
            startXYZ[2] = coordInfo.m_tabSpaceInfo.m_xyz[2];
            oneDimShape->setTabIndex(coordInfo.m_tabSpaceInfo.m_index);
            validCoordsFlag = true;
        }
        else if (windowFlag
                 && (coordInfo.m_windowSpaceInfo.m_index >= 0)) {
            startXYZ[0] = coordInfo.m_windowSpaceInfo.m_xyz[0];
            startXYZ[1] = coordInfo.m_windowSpaceInfo.m_xyz[1];
            startXYZ[2] = coordInfo.m_windowSpaceInfo.m_xyz[2];
            oneDimShape->setWindowIndex(coordInfo.m_windowSpaceInfo.m_index);
            validCoordsFlag = true;
        }
        
        if (validCoordsFlag) {
            endXYZ[0] = startXYZ[0] + diffXYZ[0];
            endXYZ[1] = startXYZ[1] + diffXYZ[1];
            endXYZ[2] = startXYZ[2] + diffXYZ[2];
            
            /*
             * Tab/Window coordinates are percentage ranging [0.0, 100.0]
             * Need to "clip" lines if they exceed the viewport's edges
             */
            const float minCoord = 1.0;
            const float maxCoord = 99.0;
            
            if (endXYZ[0] < minCoord) {
                if (diffXYZ[0] != 0.0) {
                    const float xDist = minCoord - startXYZ[0];
                    const float scaledDistance = std::fabs(xDist / diffXYZ[0]);
                    endXYZ[0] = minCoord;
                    endXYZ[1] = startXYZ[1] + (scaledDistance * diffXYZ[1]);
                }
            }
            else if (endXYZ[0] >= maxCoord) {
                if (diffXYZ[0] != 0.0) {
                    const float xDist = maxCoord - startXYZ[0];
                    const float scaledDistance = std::fabs(xDist / diffXYZ[0]);
                    endXYZ[0] = maxCoord;
                    endXYZ[1] = startXYZ[1] + (scaledDistance * diffXYZ[1]);
                }
            }
            
            if (endXYZ[1] < minCoord) {
                if (diffXYZ[1] != 0.0) {
                    const float yDist = minCoord - startXYZ[1];
                    const float scaledDistance = std::fabs(yDist / diffXYZ[1]);
                    endXYZ[1] = minCoord;
                    endXYZ[0] = startXYZ[0] + (scaledDistance * diffXYZ[0]);
                }
            }
            else if (endXYZ[1] > maxCoord) {
                if (diffXYZ[1] != 0.0) {
                    const float yDist = maxCoord - startXYZ[1];
                    const float scaledDistance = std::fabs(yDist / diffXYZ[1]);
                    endXYZ[1] = maxCoord;
                    endXYZ[0] = startXYZ[0] + (scaledDistance * diffXYZ[0]);
                }
            }
            
            oneDimShape->getStartCoordinate()->setXYZ(startXYZ);
            oneDimShape->getEndCoordinate()->setXYZ(endXYZ);
        }
    }
    
    return validCoordsFlag;
}

/**
 * Paste a multi-coordinate shape (poly-line) The first coordinate is pasted at
 * the coordinate in 'coordInfo'.
 *
 * @param multiCoordShape
 *     multi-coordinate shape that will be pasted.
 * @param coordInfo
 *     Coordinate information that will be used for the shape's 'start' coordinate.
 * @return
 *     True if the shape's coordinate was updated for pasting, else false.
 */
bool
AnnotationPasteDialog::pasteMultiCoordinateShape(AnnotationMultiCoordinateShape* multiCoordShape,
                                                 AnnotationCoordinateInformation& coordInfo)
{
    const int32_t numCoords = multiCoordShape->getNumberOfCoordinates();
    if (numCoords < 2) {
        return false;
    }
    
    bool tabFlag = false;
    bool windowFlag = false;
    
    switch (multiCoordShape->getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            tabFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            windowFlag = true;
            break;
    }
    
    bool validCoordsFlag = false;
    
    if (tabFlag
        || windowFlag) {
        float startXYZ[3];
        multiCoordShape->getCoordinate(0)->getXYZ(startXYZ);
        
        const float originalStartXYZ[3] = {
            startXYZ[0],
            startXYZ[1],
            startXYZ[2]
        };
        
        if (tabFlag
            && (coordInfo.m_tabSpaceInfo.m_index >= 0)) {
            startXYZ[0] = coordInfo.m_tabSpaceInfo.m_xyz[0];
            startXYZ[1] = coordInfo.m_tabSpaceInfo.m_xyz[1];
            startXYZ[2] = coordInfo.m_tabSpaceInfo.m_xyz[2];
            multiCoordShape->setTabIndex(coordInfo.m_tabSpaceInfo.m_index);
            validCoordsFlag = true;
        }
        else if (windowFlag
                 && (coordInfo.m_windowSpaceInfo.m_index >= 0)) {
            startXYZ[0] = coordInfo.m_windowSpaceInfo.m_xyz[0];
            startXYZ[1] = coordInfo.m_windowSpaceInfo.m_xyz[1];
            startXYZ[2] = coordInfo.m_windowSpaceInfo.m_xyz[2];
            multiCoordShape->setWindowIndex(coordInfo.m_windowSpaceInfo.m_index);
            validCoordsFlag = true;
        }
        
        if (validCoordsFlag) {
            /*
             * Set first coordinate
             */
            multiCoordShape->getCoordinate(0)->setXYZ(startXYZ);
            
            /*
             * Set additional coordinates
             */
            for (int32_t iCoord = 1; iCoord < numCoords; iCoord++) {
                float ptXYZ[3];
                multiCoordShape->getCoordinate(iCoord)->getXYZ(ptXYZ);
                
                const float diffXYZ[3] = {
                    startXYZ[0] - originalStartXYZ[0],
                    startXYZ[1] - originalStartXYZ[1],
                    startXYZ[2] - originalStartXYZ[2]
                };
                if (validCoordsFlag) {
                    ptXYZ[0] += diffXYZ[0]; // + originalStartXYZ[0];
                    ptXYZ[1] += diffXYZ[1]; //+ originalStartXYZ[1];
                    ptXYZ[2] += diffXYZ[2]; // + originalStartXYZ[2];
                    
                    /*
                     * Tab/Window coordinates are percentage ranging [0.0, 100.0]
                     * Need to "clip" lines if they exceed the viewport's edges
                     */
                    const float minCoord = 1.0;
                    const float maxCoord = 99.0;
                    
                    ptXYZ[0] = MathFunctions::limitRange(ptXYZ[0], minCoord, maxCoord);
                    ptXYZ[1] = MathFunctions::limitRange(ptXYZ[1], minCoord, maxCoord);
                    
                    multiCoordShape->getCoordinate(iCoord)->setXYZ(ptXYZ);
                }
            }
        }
    }
    
    return validCoordsFlag;
}

/**
 * Gets called when the OK button is clicked.
 */
void
AnnotationPasteDialog::okButtonClicked()
{
    AString errorMessage;
    
    bool valid = false;
    m_coordinateSelectionWidget->getSelectedCoordinateSpace(valid);
    if ( ! valid) {
        const QString msg("A coordinate space has not been selected.");
        WuQMessageBox::errorOk(this,
                               msg);
        return;
    }
    
    CaretPointer<Annotation> newAnnotation(m_annotation->clone());
    const AnnotationCoordinateSpaceEnum::Enum previousSpace = m_annotation->getCoordinateSpace();
    
    if ( ! m_coordinateSelectionWidget->setCoordinateForNewAnnotation(newAnnotation,
                                                                      errorMessage)) {
        WuQMessageBox::errorOk(this, errorMessage);
        return;
    }
    
    /*
     * Need to release annotation from its CaretPointer since the
     * annotation file will take ownership of the annotation.
     */
    Annotation* annotationPointer = newAnnotation.releasePointer();
    
    adjustTextAnnotationFontHeight(previousSpace,
                                   annotationPointer);
    
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    
    AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
    undoCommand->setModePasteAnnotation(m_annotationFile,
                                        annotationPointer);
    switch (annotationPointer->getType()) {
        case AnnotationTypeEnum::BOX:
            break;
        case AnnotationTypeEnum::BROWSER_TAB:
            CaretAssert(0);
            break;
        case AnnotationTypeEnum::COLOR_BAR:
            break;
        case AnnotationTypeEnum::IMAGE:
            break;
        case AnnotationTypeEnum::LINE:
            break;
        case AnnotationTypeEnum::OVAL:
            break;
        case AnnotationTypeEnum::POLY_LINE:
            break;
        case AnnotationTypeEnum::SCALE_BAR:
            break;
        case AnnotationTypeEnum::TEXT:
            break;
    }
    if ( ! annotationManager->applyCommand(UserInputModeEnum::Enum::ANNOTATIONS,
                                           undoCommand,
                                           errorMessage)) {
        WuQMessageBox::errorOk(this,
                               errorMessage);
    }
    
    annotationManager->selectAnnotationForEditing(m_mouseEvent.getBrowserWindowIndex(),
                                        AnnotationManager::SELECTION_MODE_SINGLE,
                                        false,
                                        annotationPointer);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendEvent(EventUserInterfaceUpdate().getPointer());
    
    m_annotationThatWasCreated = annotationPointer;

    WuQDialog::okButtonClicked();
}

/**
 * If the given annotation is a text annotation and the space is changed to 
 * surface, we may need to adjust the height if in surface montage.  This 
 * results in the surface space pixel height being the same as the tab 
 * space text pixel height when in surface montage.
 *
 * Inverse logic is need when converting text annotation from surface
 * space to another space.
 *
 * @param previousSpace
 *      Space of annotation that was on the clipboard.
 * @param annotation
 *      Annotation that is being pasted.
 */
void
AnnotationPasteDialog::adjustTextAnnotationFontHeight(const AnnotationCoordinateSpaceEnum::Enum previousSpace,
                                                      Annotation* annotation)
{
    CaretAssert(annotation);
    
    BrainOpenGLViewportContent* vpContent = m_mouseEvent.getViewportContent();
    CaretAssert(vpContent);
    BrowserTabContent* btc = vpContent->getBrowserTabContent();
    CaretAssert(btc);
    
    int32_t surfaceMontageRowCount = 1;
    const ModelSurfaceMontage* msm = btc->getDisplayedSurfaceMontageModel();
    if (msm != NULL) {
        int32_t columnCount = 1;
        msm->getSurfaceMontageNumberOfRowsAndColumns(btc->getTabNumber(),
                                                     surfaceMontageRowCount,
                                                     columnCount);
    }
    
    const AnnotationCoordinateSpaceEnum::Enum newSpace = annotation->getCoordinateSpace();
    
    const bool previousSpaceStereoOrSurfaceFlag = ((previousSpace == AnnotationCoordinateSpaceEnum::STEREOTAXIC)
                                                   || (previousSpace == AnnotationCoordinateSpaceEnum::SURFACE));
    const bool newSpaceStereoOrSurfaceFlag = ((newSpace == AnnotationCoordinateSpaceEnum::STEREOTAXIC)
                                                   || (newSpace == AnnotationCoordinateSpaceEnum::SURFACE));
    if (surfaceMontageRowCount > 1) {
        if (annotation->getType() == AnnotationTypeEnum::TEXT) {
            
            float heightMultiplier = 0.0;
            
            if ( ! previousSpaceStereoOrSurfaceFlag) {
                if (newSpaceStereoOrSurfaceFlag) {
                    /*
                     * Converting to surface
                     */
                    heightMultiplier = surfaceMontageRowCount;
                }
            }
            else {
                if ( ! newSpaceStereoOrSurfaceFlag) {
                    /*
                     * Converting from surface
                     */
                    heightMultiplier = 1.0 / surfaceMontageRowCount;
                }
            }
            
            if (heightMultiplier != 0.0) {
                AnnotationPercentSizeText* textAnn = dynamic_cast<AnnotationPercentSizeText*>(annotation);
                if (textAnn != NULL) {
                    float percentHeight = textAnn->getFontPercentViewportSize();
                    percentHeight *= heightMultiplier;
                    textAnn->setFontPercentViewportSize(percentHeight);
                }
            }
        }
    }
}


