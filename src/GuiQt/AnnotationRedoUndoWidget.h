#ifndef __ANNOTATION_REDO_UNDO_WIDGET_H__
#define __ANNOTATION_REDO_UNDO_WIDGET_H__

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


#include <stdint.h>

#include <QWidget>

#include "UserInputModeEnum.h"

namespace caret {
    class Annotation;
    class UserInputModeAnnotations;

    class AnnotationRedoUndoWidget : public QWidget {
        
        Q_OBJECT

    public:
        AnnotationRedoUndoWidget(const Qt::Orientation orientation,
                                 UserInputModeAnnotations* userInputModeAnnotations,
                                 const int32_t browserWindowIndex,
                                 QWidget* parent = 0);
        
        virtual ~AnnotationRedoUndoWidget();
        
        void updateContent(const std::vector<Annotation*>& annotations);
        
    private slots:
        void redoActionTriggered();
        
        void undoActionTriggered();
        

        // ADD_NEW_METHODS_HERE

    private:
        AnnotationRedoUndoWidget(const AnnotationRedoUndoWidget&);

        AnnotationRedoUndoWidget& operator=(const AnnotationRedoUndoWidget&);
        
        UserInputModeAnnotations* m_userInputModeAnnotations;
        
        const UserInputModeEnum::Enum m_userInputMode;
        
        const int32_t m_browserWindowIndex;
        
        QAction* m_redoAction = NULL;
        
        QAction* m_undoAction = NULL;
        
        std::vector<Annotation*> m_selectedAnnotations;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_REDO_UNDO_WIDGET_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __ANNOTATION_REDO_UNDO_WIDGET_DECLARE__

} // namespace
#endif  //__ANNOTATION_REDO_UNDO_WIDGET_H__
