#ifndef __ANNOTATION_SAMPLES_MODIFY_WIDGET_H__
#define __ANNOTATION_SAMPLES_MODIFY_WIDGET_H__

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

class QToolButton;

namespace caret {
    class Annotation;
    class AnnotationPolyhedron;

    class AnnotationSamplesModifyWidget : public QWidget {
        
        Q_OBJECT

    public:
        AnnotationSamplesModifyWidget(const Qt::Orientation orientation,
                                 const UserInputModeEnum::Enum userInputMode,
                                 const int32_t browserWindowIndex,
                                 QWidget* parent = 0);
        
        virtual ~AnnotationSamplesModifyWidget();
        
        void updateContent(const std::vector<Annotation*>& annotations);
        
    private slots:
        void lockActionToggled(bool checked);
        
        void moreActionTriggered();
        
        // ADD_NEW_METHODS_HERE

    private:
        AnnotationSamplesModifyWidget(const AnnotationSamplesModifyWidget&);

        AnnotationSamplesModifyWidget& operator=(const AnnotationSamplesModifyWidget&);
        
        const UserInputModeEnum::Enum m_userInputMode;
        
        const int32_t m_browserWindowIndex;
        
        QAction* m_lockAction = NULL;
        
        QToolButton* m_moreToolButton = NULL;
        
        QAction* m_moreAction = NULL;
        
        AnnotationPolyhedron* m_polyhedronSelected = NULL;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_SAMPLES_MODIFY_WIDGET_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __ANNOTATION_SAMPLES_MODIFY_WIDGET_DECLARE__

} // namespace
#endif  //__ANNOTATION_SAMPLES_MODIFY_WIDGET_H__
