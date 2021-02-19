#ifndef __ANNOTATION_SPATIAL_MODIFICATION_H__
#define __ANNOTATION_SPATIAL_MODIFICATION_H__

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


#include "AnnotationSizingHandleTypeEnum.h"
#include "CaretObject.h"
#include "StructureEnum.h"



namespace caret {

    class AnnotationSpatialModification : public CaretObject {
        
    public:
        AnnotationSpatialModification(const AnnotationSizingHandleTypeEnum::Enum sizingHandleType,
                                      const float viewportWidth,
                                      const float viewportHeight,
                                      const float mousePressX,
                                      const float mousePressY,
                                      const float mouseX,
                                      const float mouseY,
                                      const float mouseDX,
                                      const float mouseDY,
                                      const int32_t polyLineCoordinateIndex,
                                      const bool  startOfDraggingFlag);
        
        void setSurfaceCoordinateAtMouseXY(const StructureEnum::Enum structure,
                                     const int32_t surfaceNumberOfNodes,
                                     const int32_t surfaceNodeIndex);
        
        void setStereotaxicCoordinateAtMouseXY(const float stereotaxicX,
                                               const float stereotaxicY,
                                               const float stereotaxicZ);
        
        void setChartCoordinateAtMouseXY(const float chartX,
                                         const float chartY,
                                         const float chartZ);
        
        void setChartCoordinateAtPreviousMouseXY(const float chartX,
                                                 const float chartY,
                                                 const float chartZ);
        
        virtual ~AnnotationSpatialModification();
        

        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
    private:
        class SurfaceCoord {
        public:
            SurfaceCoord() {
                m_surfaceStructure     = StructureEnum::INVALID;
                m_surfaceNumberOfNodes = -1;
                m_surfaceNodeIndex     = -1;
                m_surfaceNodeValid     = false;
            }
            StructureEnum::Enum m_surfaceStructure;
            
            int32_t m_surfaceNumberOfNodes;
            
            int32_t m_surfaceNodeIndex;
            
            bool m_surfaceNodeValid;
        };
        
        class ChartCoord {
        public:
            ChartCoord() {
                m_chartXYZ[0] = 0.0;
                m_chartXYZ[1] = 0.0;
                m_chartXYZ[2] = 0.0;
                
                m_chartXYZValid = false;
            }
            
            float m_chartXYZ[3];
            
            bool m_chartXYZValid;
        };
        
        
        class StereotaxicCoord {
        public:
            StereotaxicCoord() {
                m_stereotaxicXYZ[0] = 0.0;
                m_stereotaxicXYZ[1] = 0.0;
                m_stereotaxicXYZ[2] = 0.0;
                
                m_stereotaxicValid = false;
            }
            
            float m_stereotaxicXYZ[3];
            
            bool m_stereotaxicValid;
        };
        
        AnnotationSpatialModification(const AnnotationSpatialModification&);

        AnnotationSpatialModification& operator=(const AnnotationSpatialModification&);
        
        const AnnotationSizingHandleTypeEnum::Enum m_sizingHandleType;
        
        const float m_viewportWidth;
        
        const float m_viewportHeight;
        
        const float m_mousePressX;
        
        const float m_mousePressY;
        
        const float m_mouseX;
        
        const float m_mouseY;
        
        const float m_mouseDX;
        
        const float m_mouseDY;
        
        const int32_t m_polyLineCoordinateIndex;
        
        const bool  m_startOfDraggingFlag;
        
        ChartCoord m_chartCoordAtMouseXY;
        
        ChartCoord m_chartCoordAtPreviousMouseXY;
        
        SurfaceCoord m_surfaceCoordinateAtMouseXY;
        
        StereotaxicCoord m_stereotaxicCoordinateAtMouseXY;
        
        // ADD_NEW_MEMBERS_HERE

        friend class AnnotationMultiCoordinateShape;
        friend class AnnotationText;
        friend class AnnotationOneCoordinateShape;
        friend class AnnotationTwoCoordinateShape;
    };
    
#ifdef __ANNOTATION_SPATIAL_MODIFICATION_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __ANNOTATION_SPATIAL_MODIFICATION_DECLARE__

} // namespace
#endif  //__ANNOTATION_SPATIAL_MODIFICATION_H__
