#ifndef __PALETTE_PIXMAP_PAINTER_H__
#define __PALETTE_PIXMAP_PAINTER_H__

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



#include <memory>

#include <QPixmap>

#include "CaretObject.h"

namespace caret {

    class Palette;
    class PaletteNew;
    
    class PalettePixmapPainter : public CaretObject {
        
    public:
        enum class Mode {
            INTERPOLATE_OFF,
            INTERPOLATE_ON,
            INTERPOLATE_ON_LINES_AT_SCALARS,
            SCALAR_LINES
        };
        
        PalettePixmapPainter(const Palette* palette,
                             const Mode mode);
        
        PalettePixmapPainter(const Palette* palette,
                             const QSize& pixmapSize,
                             const Mode mode);
        
        PalettePixmapPainter(const PaletteNew* palette,
                             const QSize& pixmapSize,
                             const Mode mode);
        
        virtual ~PalettePixmapPainter();
        
        PalettePixmapPainter(const PalettePixmapPainter&) = delete;

        PalettePixmapPainter& operator=(const PalettePixmapPainter&) = delete;

        QPixmap getPixmap() const;

        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
    private:
        enum class DrawScalarLinesMode {
            LINES_OFF,
            LINES_ON
        };
        
        void createPalettePixmapInterpolateOff(const Palette* palette,
                                               const qreal pixmapWidth,
                                               const qreal pixmapHeight);
        
        void createPalettePixmapInterpolateOn(const Palette* palette,
                                              const qreal pixmapWidth,
                                              const qreal pixmapHeight);
        
        float paletteNonInterpScalarToPixMap(const float pixmapWidth,
                                             const float scalarIn);

        void createPalettePixmapInterpolateOff(const PaletteNew* palette,
                                               const qreal pixmapWidth,
                                               const qreal pixmapHeight);
        
        void createPalettePixmapInterpolateOn(const PaletteNew* palette,
                                              const qreal pixmapWidth,
                                              const qreal pixmapHeight,
                                              const DrawScalarLinesMode drawLinesMode);
        
        void createPalettePixmapScalarLines(const PaletteNew* palette,
                                            const qreal pixmapWidth,
                                            const qreal pixmapHeight);
        
        void drawScalarLines(const PaletteNew* palette,
                             QPainter& painter,
                             QPen& pen,
                             QPixmap& pixmap,
                             const qreal lineWidth);

        void drawLineInColorBar(QPainter& painter,
                                QPen& pen,
                                const float rgbFloat[3],
                                const qreal x,
                                const qreal y,
                                const qreal pixmapHeight,
                                const int32_t lineWidth);

        const Mode m_mode = Mode::INTERPOLATE_ON;
        
        QPixmap m_pixmap;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __PALETTE_PIXMAP_PAINTER_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __PALETTE_PIXMAP_PAINTER_DECLARE__

} // namespace
#endif  //__PALETTE_PIXMAP_PAINTER_H__
