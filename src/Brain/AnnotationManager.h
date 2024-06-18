#ifndef __ANNOTATION_MANAGER_H__
#define __ANNOTATION_MANAGER_H__

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

#include <memory>

#include "AnnotationAlignmentEnum.h"
#include "AnnotationAndFile.h"
#include "AnnotationCoordinateSpaceEnum.h"
#include "AnnotationDistributeEnum.h"
#include "AnnotationGroupingModeEnum.h"
#include "AnnotationStackingOrderTypeEnum.h"
#include "BrainConstants.h"
#include "CaretObject.h"
#include "CaretUndoCommand.h"
#include "CaretPointer.h"
#include "EventListenerInterface.h"
#include "SceneableInterface.h"
#include "UserInputModeEnum.h"

namespace caret {
    class Annotation;
    class AnnotationArrangerInputs;
    class AnnotationBrowserTab;
    class AnnotationClipboard;
    class AnnotationFile;
    class AnnotationGroupKey;
    class AnnotationPolyhedron;
    class AnnotationRedoUndoCommand;
    class AnnotationEditingSelectionInformation;
    class Brain;
    class BrowserTabContent;
    class CaretUndoStack;
    class EventGetDisplayedDataFiles;
    class SceneClassAssistant;

    class AnnotationManager : public CaretObject, public EventListenerInterface, public SceneableInterface {
        
    public:
        /**
         * Selection mode with the comments and functionality similar to QAbstractItemView::SelectionMode
         */
        enum SelectionMode {
            /*
             * When the user selects an annotation, any already-selected annotation becomes unselected,
             * and the user cannot unselect the selected annotation by clicking on it.
             */
            SELECTION_MODE_SINGLE,
            /**
             * When the user selects an annotation in the usual way, the selection is cleared
             * and the new annotation selected. However, if the user presses the Ctrl key when
             * clicking on an annotation, the clicked annotation gets toggled and all other annotation
             * are left untouched.
             */
            SELECTION_MODE_EXTENDED
        };
        
        AnnotationManager(const UserInputModeEnum::Enum userInputMode,
                          Brain* brain);
        
        virtual ~AnnotationManager();

        bool applyCommand(AnnotationRedoUndoCommand* command,
                          AString& errorMessageOut);
        
        bool applyCommandInWindow(AnnotationRedoUndoCommand* command,
                                  const int32_t windowIndex,
                                  AString& errorMessageOut);
        
        void reset();
        
        void deselectAllAnnotationsForEditing(const int32_t windowIndex);
        
        void selectAnnotationForEditing(const int32_t windowIndex,
                              const SelectionMode selectionMode,
                              const bool shiftKeyDownFlag,
                              Annotation* selectedAnnotation);
        
        void setAnnotationsForEditing(const int32_t windowIndex,
                                      const std::vector<Annotation*>& selectedAnnotations);
        
        bool isAnnotationSelectedForEditingDeletable(const int32_t windowIndex) const;
        
        std::vector<Annotation*> getAllAnnotations() const;
        
        std::vector<Annotation*> getAnnotationsDrawnInSameWindowAndSpace(const Annotation* annotation,
                                                                         const int32_t windowIndex) const;
        
        const AnnotationEditingSelectionInformation* getAnnotationEditingSelectionInformation(const int32_t windowIndex) const;
        
        std::vector<Annotation*> getAnnotationsSelectedForEditing(const int32_t windowIndex) const;
        
        std::vector<Annotation*> getAnnotationsSelectedForEditingInSpaces(const int32_t windowIndex,
                                                                const std::vector<AnnotationCoordinateSpaceEnum::Enum>& spaces) const;
        
        void getAnnotationsAndFilesSelectedForEditing(const int32_t windowIndex,
                                                      std::vector<AnnotationAndFile>& annotationsAndFileOut) const;
        
        void getAnnotationsAndFilesSelectedForEditingIncludingLabels(const int32_t windowIndex,
                                                                     std::vector<AnnotationAndFile>& annotationsAndFileOut) const;
        
        AnnotationClipboard* getClipboard();
        
        const AnnotationClipboard* getClipboard() const;
        
        CaretUndoStack* getCommandRedoUndoStack();
        
        void getDisplayedAnnotationFiles(EventGetDisplayedDataFiles* displayedFilesEvent,
                                         std::vector<AnnotationFile*>& displayedAnnotationFilesOut) const;
        
        bool alignAnnotations(const AnnotationArrangerInputs& arrangerInputs,
                              const AnnotationAlignmentEnum::Enum alignment,
                              AString& errorMessageOut);
        
        bool distributeAnnotations(const AnnotationArrangerInputs& arrangerInputs,
                                   const AnnotationDistributeEnum::Enum distribute,
                                   AString& errorMessageOut);

        bool applyGroupingMode(const int32_t windowIndex,
                               const AnnotationGroupingModeEnum::Enum groupingMode,
                               AString& errorMessageOut);
        
        bool isGroupingModeValid(const int32_t windowIndex,
                                 const AnnotationGroupingModeEnum::Enum groupingMode) const;
        
        bool applyStackingOrder(const std::vector<Annotation*>& annotations,
                                const Annotation* selectedAnnotation,
                                const AnnotationStackingOrderTypeEnum::Enum orderType,
                                const int32_t windowIndex,
                                AString& errorMessageOut);
        
        bool moveTabOrWindowAnnotationToFront(Annotation* annotation,
                                              AString& errorMessageOut);
        
        bool shrinkAndExpandSelectedBrowserTabAnnotation(const std::vector<BrowserTabContent*>& tabsInWindow,
                                                         const int32_t windowIndex,
                                                         AString& errorMessageOut);

        std::vector<Annotation*> getAnnotationsInSameSpace(const Annotation* annotation);
        
        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
        virtual void receiveEvent(Event* event);

        virtual SceneClass* saveToScene(const SceneAttributes* sceneAttributes,
                                        const AString& instanceName);

        virtual void restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass);
          
// If there will be sub-classes of this class that need to save
// and restore data from scenes, these pure virtual methods can
// be uncommented to force their implemetation by sub-classes.
//    protected: 
//        virtual void saveSubClassDataToScene(const SceneAttributes* sceneAttributes,
//                                             SceneClass* sceneClass) = 0;
//
//        virtual void restoreSubClassDataFromScene(const SceneAttributes* sceneAttributes,
//                                                  const SceneClass* sceneClass) = 0;

    private:
        AnnotationManager(const AnnotationManager&);

        AnnotationManager& operator=(const AnnotationManager&);
        
        void processExtendedModeSelectionForEditing(const int32_t windowIndex,
                                          const bool shiftKeyDownFlag,
                                          Annotation* selectedAnnotation);
        
        void processSingleModeSelectionForEditing(const int32_t windowIndex,
                                        Annotation* selectedAnnotation);
        
        const UserInputModeEnum::Enum m_userInputMode;
        
        /** Brain owning this manager */
        Brain* m_brain;
        
        SceneClassAssistant* m_sceneAssistant;
        
        /*
         * DO NOT directly reference this.  Instead, call this class'
         * getAnnotationEditingSelectionInformation() method so that the selection 
         * information is updated.
         */
        mutable AnnotationEditingSelectionInformation* m_selectionInformation[BrainConstants::MAXIMUM_NUMBER_OF_BROWSER_WINDOWS];
        
        std::unique_ptr<AnnotationClipboard> m_clipboard;
        
        CaretPointer<CaretUndoStack> m_annotationsExceptBrowserTabsRedoUndoStack;
        
        CaretPointer<CaretUndoStack> m_browserTabAnnotationsRedoUndoStack;
        
        CaretPointer<CaretUndoStack> m_samplesAnnotationsRedoUndoStack;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_MANAGER_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __ANNOTATION_MANAGER_DECLARE__

} // namespace
#endif  //__ANNOTATION_MANAGER_H__
