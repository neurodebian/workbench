#ifndef __WB_MACRO_CUSTOM_OPERATION_SURFACE_DEFAULT_COLOR_H__
#define __WB_MACRO_CUSTOM_OPERATION_SURFACE_DEFAULT_COLOR_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2019 Washington University School of Medicine
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

#include "WbMacroCustomOperationBase.h"



namespace caret {

    class BrowserTabContent;
    
    class WbMacroCustomOperationSurfaceDefaultColor : public WbMacroCustomOperationBase {
        
    public:
        WbMacroCustomOperationSurfaceDefaultColor();
        
        virtual ~WbMacroCustomOperationSurfaceDefaultColor();
        
        WbMacroCustomOperationSurfaceDefaultColor(const WbMacroCustomOperationSurfaceDefaultColor&) = delete;

        WbMacroCustomOperationSurfaceDefaultColor& operator=(const WbMacroCustomOperationSurfaceDefaultColor&) = delete;
                
        virtual bool executeCommand(QWidget* parent,
                                    const WuQMacroExecutorMonitor* executorMonitor,
                                    const WuQMacroExecutorOptions* executorOptions,
                                    const WuQMacroCommand* macroCommand) override;
        
        virtual WuQMacroCommand* createCommand() override;
        

        // ADD_NEW_METHODS_HERE

    private:
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __WB_MACRO_CUSTOM_OPERATION_SURFACE_DEFAULT_COLOR_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __WB_MACRO_CUSTOM_OPERATION_SURFACE_DEFAULT_COLOR_DECLARE__

} // namespace
#endif  //__WB_MACRO_CUSTOM_OPERATION_SURFACE_DEFAULT_COLOR_H__
