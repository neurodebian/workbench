#ifndef __ANNOTATION_COORDINATE_SELECTION_WIDGET_H__
#define __ANNOTATION_COORDINATE_SELECTION_WIDGET_H__

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
#include <QWidget>

//#include "AnnotationCoordinateInformation.h"
#include "AnnotationCoordinateSpaceEnum.h"
#include "AnnotationTypeEnum.h"
#include "UserInputModeAnnotations.h"

class QButtonGroup;
class QRadioButton;

namespace caret {

    class Annotation;
    class AnnotationCoordinate;
    class AnnotationImage;
    
    class AnnotationCoordinateSelectionWidget : public QWidget {
        
        Q_OBJECT

    public:
        AnnotationCoordinateSelectionWidget(const AnnotationTypeEnum::Enum annotationType,
                                            const AnnotationCoordinateInformation& coordInfo,
                                            const AnnotationCoordinateInformation* optionalSecondCoordInfo,
                                            QWidget* parent = 0);
        
        virtual ~AnnotationCoordinateSelectionWidget();
        
        void selectCoordinateSpace(const AnnotationCoordinateSpaceEnum::Enum coordSpace);

        AnnotationCoordinateSpaceEnum::Enum getSelectedCoordinateSpace(bool& validOut) const;
        
        bool changeAnnotationCoordinate(Annotation* annotation,
                                        QString& errorMessageOut);
        
        bool setCoordinateForNewAnnotation(Annotation* annotation,
                                           QString& errorMessageOut);
        
        // ADD_NEW_METHODS_HERE

    private:
        AnnotationCoordinateSelectionWidget(const AnnotationCoordinateSelectionWidget&);

        AnnotationCoordinateSelectionWidget& operator=(const AnnotationCoordinateSelectionWidget&);
        
        void updateAnnotationDisplayProperties(const Annotation* annotation);
        
        void setWidthAndHeightForImage(AnnotationImage* imageAnn);
        
        QRadioButton* createRadioButtonForSpace(const AnnotationCoordinateSpaceEnum::Enum space);

        QButtonGroup* m_spaceButtonGroup;
        
        const AnnotationTypeEnum::Enum m_annotationType;
        
        const AnnotationCoordinateInformation& m_coordInfo;
        
        const AnnotationCoordinateInformation* m_optionalSecondCoordInfo;
        
        std::vector<std::unique_ptr<AnnotationCoordinateInformation>> m_optionalMultiCoordInfo;
        
        static const QString s_SPACE_PROPERTY_NAME;
        

        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __ANNOTATION_COORDINATE_SELECTION_WIDGET_DECLARE__
    const QString AnnotationCoordinateSelectionWidget::s_SPACE_PROPERTY_NAME = "SPACE_NAME";
#endif // __ANNOTATION_COORDINATE_SELECTION_WIDGET_DECLARE__

} // namespace
#endif  //__ANNOTATION_COORDINATE_SELECTION_WIDGET_H__
