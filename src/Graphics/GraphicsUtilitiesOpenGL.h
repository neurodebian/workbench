#ifndef __GRAPHICS_UTILITIES_OPEN_G_L_H__
#define __GRAPHICS_UTILITIES_OPEN_G_L_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2017 Washington University School of Medicine
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
#include <stack>

#include "CaretObject.h"
#include "CaretOpenGLInclude.h"

namespace caret {

    class GraphicsOpenGLError;
    
    class GraphicsUtilitiesOpenGL : public CaretObject {
        
    public:
        
        static float convertMillimetersToPercentageOfViewportHeight(const float millimeters);
        
        static float convertPercentageOfViewportHeightToMillimeters(const float percentOfViewportHeight);
        
        static float convertPercentageOfViewportHeightToPixels(const float percentOfViewportHeight);
        
        static float convertMillimetersToPixels(const float millimeters);
        
        static float convertPixelsToPercentageOfViewportHeight(const float pixels);
        
        static float convertPixelsToMillimeters(const float pixels);
        
        static void getTextCompressionEnumInfo(const GLenum enumValue,
                                               QString& nameOut,
                                               QString& decimalValueOut,
                                               QString& hexadecimalValueOut);
        
        static void resetOpenGLError();
        
        static std::unique_ptr<GraphicsOpenGLError> getOpenGLError(const AString& message = "");
        
        static void printMatrixDepths(const AString messagePrefix);
        
        static bool isVersionOrGreater(const int32_t majorVersion,
                                       const int32_t minorVersion);
        
        static int32_t getMajorVersion();
        
        static int32_t getMinorVersion();
        
        static QString getVersion();
        
        static int32_t getTextureWidthHeightMaximumDimension();
        
        static int32_t getTextureDepthMaximumDimension();
        
        static void pushMatrix();
        
        static void popMatrix();

        static bool unproject(const float windowX,
                              const float windowY,
                              float modelXyzOut[3]);
        
    private:
        GraphicsUtilitiesOpenGL();
        
        virtual ~GraphicsUtilitiesOpenGL();
        
        GraphicsUtilitiesOpenGL(const GraphicsUtilitiesOpenGL&);

        GraphicsUtilitiesOpenGL& operator=(const GraphicsUtilitiesOpenGL&);

        static void setMajorMinorVersion(const int32_t majorVersion,
                                         const int32_t minorVersion);
        
        static void setMaximumTextureDimension(const int32_t widthHeightMaximumDimension,
                                               const int32_t depthMaximumDimension);
        
        static int32_t s_majorVersion;
        
        static int32_t s_minorVersion;
        
        static int32_t s_textureWidthHeightMaximumDimension;
        
        static int32_t s_textureDepthMaximumDimension;
        
        static std::stack<std::array<double, 16>> s_projectionMatrixStack;
        
        friend class BrainOpenGL;
    };
    
#ifdef __GRAPHICS_UTILITIES_OPEN_G_L_DECLARE__
    int32_t GraphicsUtilitiesOpenGL::s_majorVersion = 0;
    int32_t GraphicsUtilitiesOpenGL::s_minorVersion = 0;
    
    /*
     * Default texture dimensions to typical values
     */
    int32_t GraphicsUtilitiesOpenGL::s_textureWidthHeightMaximumDimension = 16384;
    int32_t GraphicsUtilitiesOpenGL::s_textureDepthMaximumDimension = 2048;
    
    std::stack<std::array<double, 16>> GraphicsUtilitiesOpenGL::s_projectionMatrixStack;
    
#endif // __GRAPHICS_UTILITIES_OPEN_G_L_DECLARE__

} // namespace
#endif  //__GRAPHICS_UTILITIES_OPEN_G_L_H__
