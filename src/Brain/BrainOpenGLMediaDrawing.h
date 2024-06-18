#ifndef __BRAIN_OPEN_G_L_MEDIA_DRAWING_H__
#define __BRAIN_OPEN_G_L_MEDIA_DRAWING_H__

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


#include <array>
#include <memory>

#include "CaretObject.h"
#include "GraphicsTextureMagnificationFilterEnum.h"
#include "GraphicsTextureMinificationFilterEnum.h"


namespace caret {

    class BrainOpenGLFixedPipeline;
    class BrainOpenGLViewportContent;
    class BrowserTabContent;
    class GraphicsObjectToWindowTransform;
    class GraphicsPrimitiveV3fT2f;
    class MediaFile;
    class MediaOverlaySet;
    class ModelMedia;
    
    class BrainOpenGLMediaDrawing : public CaretObject {
        
    public:
        static GraphicsTextureMagnificationFilterEnum::Enum getTextureMagnificationFilter();
        
        static GraphicsTextureMinificationFilterEnum::Enum  getTextureMinificationFilter();
        
        static void setTextureMagnificationFilter(const GraphicsTextureMagnificationFilterEnum::Enum magFilter);
        
        static void setTextureMinificationFilter(const GraphicsTextureMinificationFilterEnum::Enum minFilter);
        
        BrainOpenGLMediaDrawing();
        
        virtual ~BrainOpenGLMediaDrawing();
        
        void draw(BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                  const BrainOpenGLViewportContent* viewportContent,
                  BrowserTabContent* browserTabContent,
                  ModelMedia* mediaModel,
                  const std::array<int32_t, 4>& viewport);
        
        BrainOpenGLMediaDrawing(const BrainOpenGLMediaDrawing&) = delete;

        BrainOpenGLMediaDrawing& operator=(const BrainOpenGLMediaDrawing&) = delete;
        

        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
    private:
        void drawModelLayers(const BrainOpenGLViewportContent* viewportContent,
                             const GraphicsObjectToWindowTransform* transform,
                             const int32_t tabIndex,
                             const float viewportHeight);
        
        void processMediaFileSelection(const int32_t tabIndex,
                                       const int32_t overlayIndex,
                                       MediaFile* mediaFile,
                                       const int32_t selectedFrameIndex,
                                       const bool allFramesSelectedFlag,
                                       GraphicsPrimitiveV3fT2f* primitive);
        
        bool getOrthoBounds(MediaOverlaySet* mediaOverlaySet,
                            double& orthoLeftOut,
                            double& orthoRightOut,
                            double& orthoBottomOut,
                            double& orthoTopOut);

        BrainOpenGLFixedPipeline* m_fixedPipelineDrawing = NULL;
        
        BrowserTabContent* m_browserTabContent = NULL;
        
        ModelMedia* m_mediaModel;
        
        std::array<int32_t, 4> m_viewport;
        
        static GraphicsTextureMagnificationFilterEnum::Enum s_textureMagnificationFilter;
        static GraphicsTextureMinificationFilterEnum::Enum  s_textureMinificationFilter;
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __BRAIN_OPEN_G_L_MEDIA_DRAWING_DECLARE__
    GraphicsTextureMagnificationFilterEnum::Enum BrainOpenGLMediaDrawing::s_textureMagnificationFilter = GraphicsTextureMagnificationFilterEnum::LINEAR;
    GraphicsTextureMinificationFilterEnum::Enum  BrainOpenGLMediaDrawing::s_textureMinificationFilter  = GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR;
#endif // __BRAIN_OPEN_G_L_MEDIA_DRAWING_DECLARE__

} // namespace
#endif  //__BRAIN_OPEN_G_L_MEDIA_DRAWING_H__
