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

#include "AnnotationAlignmentEnum.h"
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
    class AnnotationFile;
    class AnnotationGroupKey;
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
        
        AnnotationManager(Brain* brain);
        
        virtual ~AnnotationManager();

        bool applyCommand(const UserInputModeEnum::Enum userInputMode,
                          AnnotationRedoUndoCommand* command,
                          AString& errorMessageOut);
        
        bool applyCommandInWindow(const UserInputModeEnum::Enum userInputMode,
                                  AnnotationRedoUndoCommand* command,
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
                                    std::vector<std::pair<Annotation*, AnnotationFile*> >& annotationsAndFileOut) const;
        
        void getAnnotationsAndFilesSelectedForEditingIncludingLabels(const int32_t windowIndex,
                                              std::vector<std::pair<Annotation*, AnnotationFile*> >& annotationsAndFileOut) const;
        
        bool isAnnotationOnClipboardValid() const;
        
        AnnotationFile* getAnnotationFileOnClipboard() const;
        
        const Annotation* getAnnotationOnClipboard() const;
        
        Annotation* getCopyOfAnnotationOnClipboard() const;
        
        void copyAnnotationToClipboard(const AnnotationFile* annotationFile,
                                       const Annotation* annotation);
        
        const Annotation* getAnnotationBeingDrawnInWindow(const int32_t windowIndex) const;
        
        void setAnnotationBeingDrawnInWindow(const int32_t windowIndex,
                                             const Annotation* annotation);
        
        CaretUndoStack* getCommandRedoUndoStack(const UserInputModeEnum::Enum userInputMode);
        
        void getDisplayedAnnotationFiles(EventGetDisplayedDataFiles* displayedFilesEvent,
                                         std::vector<AnnotationFile*>& displayedAnnotationFilesOut) const;
        
        bool alignAnnotations(const UserInputModeEnum::Enum userInputMode,
                              const AnnotationArrangerInputs& arrangerInputs,
                              const AnnotationAlignmentEnum::Enum alignment,
                              AString& errorMessageOut);
        
        bool distributeAnnotations(const UserInputModeEnum::Enum userInputMode,
                                   const AnnotationArrangerInputs& arrangerInputs,
                                   const AnnotationDistributeEnum::Enum distribute,
                                   AString& errorMessageOut);

        bool applyGroupingMode(const UserInputModeEnum::Enum userInputMode,
                               const int32_t windowIndex,
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
                                                         const UserInputModeEnum::Enum userInputMode,
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
        
        SceneClassAssistant* m_sceneAssistant;

        /** Brain owning this manager */
        Brain* m_brain;
        
        Annotation* m_annotationBeingDrawnInWindow[BrainConstants::MAXIMUM_NUMBER_OF_BROWSER_WINDOWS];
        
        /*
         * DO NOT directly reference this.  Instead, call this class'
         * getAnnotationEditingSelectionInformation() method so that the selection 
         * information is updated.
         */
        mutable AnnotationEditingSelectionInformation* m_selectionInformation[BrainConstants::MAXIMUM_NUMBER_OF_BROWSER_WINDOWS];
        
        /** 
         * Do not use a Caret Pointer for this as it points to a file in the brain.
         * If a pointer was used it may get deleted which will cause deletion of the
         * file that is owned by the Brain.
         */
        mutable AnnotationFile* m_clipboardAnnotationFile;
        
        CaretPointer<Annotation> m_clipboardAnnotation;
        
        CaretPointer<CaretUndoStack> m_annotationsExceptBrowserTabsRedoUndoStack;
        
        CaretPointer<CaretUndoStack> m_browserTabAnnotationsRedoUndoStack;
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_MANAGER_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __ANNOTATION_MANAGER_DECLARE__

} // namespace
#endif  //__ANNOTATION_MANAGER_H__
