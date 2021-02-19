
/*LICENSE_START*/
/*
 *  Copyright (C) 2014  Washington University School of Medicine
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

#define __CARET_COLOR_TOOLBUTTON_DECLARE__
#include "CaretColorToolButton.h"
#undef __CARET_COLOR_TOOLBUTTON_DECLARE__

#include <QComboBox>
#include <QColorDialog>
#include <QHBoxLayout>
#include <QLabel>

#include "CaretAssert.h"
#include "CaretColorEnumMenu.h"
#include "CaretLogger.h"
#include "WuQMacroManager.h"
#include "WuQtUtilities.h"

using namespace caret;


    
/**
 * \class caret::CaretColorToolButton
 * \brief Control for selection of Caret Color enumerated types.
 *
 * Control for selection of Caret Color enumerated types.
 */

/**
 * Constructor.
 *
 * @param parent
 *     Parent widget.
 */
CaretColorToolButton::CaretColorToolButton(QWidget* parent)
: CaretColorToolButton(CustomColorMode::DISABLED,
                       NoneColorMode::DISABLED,
                       parent)
{
    /* delegating constructor above */
}

/**
 * Constructor.
 *
 * @param customColorMode
 *     Mode for custom color
 * @param noneColorMode
 *     Mode for none color
 * @param parent
 *     Parent object.
 */
CaretColorToolButton::CaretColorToolButton(const CustomColorMode customColorMode,
                                           const NoneColorMode noneColorMode,
                                           QWidget* parent)
: QToolButton(parent),
m_customColorMode(customColorMode),
m_noneColorMode(noneColorMode)
{
    int64_t caretColorOptions = 0;
    switch (m_customColorMode) {
        case CustomColorMode::DISABLED:
            break;
        case CustomColorMode::EDITABLE:
            caretColorOptions |= CaretColorEnum::OPTION_INCLUDE_CUSTOM_COLOR;
            break;
        case CustomColorMode::FIXED:
            caretColorOptions |= CaretColorEnum::OPTION_INCLUDE_CUSTOM_COLOR;
            break;
    }
    
    switch (m_noneColorMode) {
        case NoneColorMode::DISABLED:
            break;
        case NoneColorMode::ENABLED:
            caretColorOptions |= CaretColorEnum::OPTION_INCLUDE_NONE_COLOR;
            break;
    }
    m_caretColorEnumMenu = new CaretColorEnumMenu(caretColorOptions);
        
    QObject::connect(this, &QToolButton::clicked,
                     this, &CaretColorToolButton::toolButtonClicked);
    
    m_caretColor.setCaretColorEnum(CaretColorEnum::BLACK);
}

/**
 * Called when the button is clicked to pop-up a menu for color selection
 */
void
CaretColorToolButton::toolButtonClicked()
{
    setSelectedColor(m_caretColor);
    QAction* action = m_caretColorEnumMenu->exec(this->mapToGlobal(QPoint(0,0)));
    if (action != NULL) {
        CaretColorEnum::Enum colorEnum = m_caretColorEnumMenu->getSelectedColor();
        caretColorMenuSelected(colorEnum);
    }
}

/**
 * Set the icon color
 * @param rgba
 * RGBA components
 */
void
CaretColorToolButton::updateIconColor()
{
    const std::array<float, 4> rgbaF = m_caretColor.getFloatRGBA();
    QPixmap pm(WuQtUtilities::createCaretColorEnumPixmap(this, 10, 10,
                                                         CaretColorEnum::CUSTOM,
                                                         rgbaF.data(), false));
    QIcon icon(pm);
    this->setIcon(icon);
}

/**
 * Destructor.
 */
CaretColorToolButton::~CaretColorToolButton()
{
    
}

/**
 * @return The selected color.
 */
CaretColor
CaretColorToolButton::getSelectedColor()
{
    return m_caretColor;
}

/**
 * Set the selected color.
 * @param color
 *   New color for selection.
 */
void 
CaretColorToolButton::setSelectedColor(const CaretColor& color)
{
    m_caretColor = color;
    m_caretColorEnumMenu->setSelectedColor(m_caretColor.getCaretColorEnum());
    m_caretColorEnumMenu->setCustomIconColor(m_caretColor.getCustomColorFloatRGBA().data());
    updateIconColor();
}

/**
 * Called when a color is selected.
 * @param colorEnum
 *   Color selected
 */
void 
CaretColorToolButton::caretColorMenuSelected(const CaretColorEnum::Enum colorEnum)
{
    if (colorEnum == CaretColorEnum::CUSTOM) {
        switch (m_customColorMode) {
            case CustomColorMode::DISABLED:
                break;
            case CustomColorMode::EDITABLE:
            {
                std::array<uint8_t, 4> rgba = m_caretColor.getCustomColorRGBA();
                const QColor initialColor(rgba[0],
                                          rgba[1],
                                          rgba[2]);
                
                QColorDialog colorDialog(this);
                colorDialog.setOption(QColorDialog::DontUseNativeDialog);
                colorDialog.setWindowTitle("Choose Color");
                colorDialog.setCurrentColor(initialColor);
                
                if (colorDialog.exec() == QColorDialog::Accepted) {
                    const QColor newColor = colorDialog.currentColor();
                    rgba[0] = newColor.red();
                    rgba[1] = newColor.green();
                    rgba[2] = newColor.blue();
                    m_caretColor.setCustomColorRGBA(rgba);
                }
                else {
                    return;
                }
            }
                break;
            case CustomColorMode::FIXED:
                break;
        }
    }
    
    m_caretColor.setCaretColorEnum(colorEnum);
    updateIconColor();
    emit colorSelected(m_caretColor);
}

