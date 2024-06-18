#ifndef __ANNOTATION_INSERT_NEW_WIDGET_H__
#define __ANNOTATION_INSERT_NEW_WIDGET_H__

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


#include <QWidget>

#include "AnnotationCoordinateSpaceEnum.h"
#include "AnnotationTypeEnum.h"
#include "EventListenerInterface.h"
#include "UserInputModeEnum.h"

class QAbstractButton;
class QActionGroup;
class QSpinBox;
class QRadioButton;
class QToolButton;

namespace caret {
    class AnnotationMenuFileSelection;
    
    class AnnotationInsertNewWidget : public QWidget, public EventListenerInterface {
        
        Q_OBJECT

    public:
        AnnotationInsertNewWidget(const UserInputModeEnum::Enum userInputMode,
                                  const int32_t browserWindowIndex,
                                  QWidget* parent = 0);
        
        virtual ~AnnotationInsertNewWidget();
        
        virtual void receiveEvent(Event* event);

        // ADD_NEW_METHODS_HERE

        void updateContent();
        
    private slots:
        void itemSelectedFromFileSelectionMenu();
        
        void spaceOrShapeActionTriggered();
        
        void newSampleActionTriggered();
        
        void newSampleDepthValueChanged(int value);
        
    private:
        enum class WidgetMode {
            INVALID,
            ANNOTATIONS,
            SAMPLES
        };
        
        AnnotationInsertNewWidget(const AnnotationInsertNewWidget&);

        AnnotationInsertNewWidget& operator=(const AnnotationInsertNewWidget&);
        
        QToolButton* createShapeToolButton(const AnnotationTypeEnum::Enum annotationType,
                                           QActionGroup* actionGroup);
        
        QToolButton* createSpaceToolButton(const AnnotationCoordinateSpaceEnum::Enum annotationSpace,
                                           QActionGroup* actionGroup);
        
        QPixmap createShapePixmap(const QWidget* widget,
                                  const AnnotationTypeEnum::Enum annotationType);
        
        QPixmap createSpacePixmap(const QWidget* widget,
                                  const AnnotationCoordinateSpaceEnum::Enum annotationSpace);
        
        QToolButton* createFileSelectionToolButton();
        
        void enableDisableSpaceActions();
        
        void enableDisableShapeActions();
        
        void createAnnotationsWidgets();
        
        void createEditSamplesWidgets();
        
        const UserInputModeEnum::Enum m_userInputMode;
        
        const int32_t m_browserWindowIndex;
        
        WidgetMode m_widgetMode = WidgetMode::INVALID;
        
        QActionGroup* m_spaceActionGroup;
        
        QActionGroup* m_shapeActionGroup;
        
        QAction* m_fileSelectionToolButtonAction;
        
        AnnotationMenuFileSelection* m_fileSelectionMenu;
        
        QToolButton* m_polygonToolButton = NULL;
        
        QToolButton* m_polyLineToolButton = NULL;
        
        QAction* m_newSampleAction = NULL;
        
        QSpinBox* m_newSampleDepthSpinBox = NULL;
        
        int32_t m_previousNewSampleDepthSpinBoxValue = 3;
        
        static AString s_previousImageFileDirectory;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_INSERT_NEW_WIDGET_DECLARE__
    AString AnnotationInsertNewWidget::s_previousImageFileDirectory;
#endif // __ANNOTATION_INSERT_NEW_WIDGET_DECLARE__

} // namespace
#endif  //__ANNOTATION_INSERT_NEW_WIDGET_H__
