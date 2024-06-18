
/*LICENSE_START*/
/*
 *  Copyright (C) 2016 Washington University School of Medicine
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

#define __DISPLAY_GROUP_AND_TAB_ITEM_VIEW_CONTROLLER_DECLARE__
#include "DisplayGroupAndTabItemViewController.h"
#undef __DISPLAY_GROUP_AND_TAB_ITEM_VIEW_CONTROLLER_DECLARE__

#include <QLabel>
#include <QMenu>
#include <QTreeWidget>
#include <QToolButton>
#include <QVBoxLayout>

#include "Annotation.h"
#include "AnnotationGroup.h"
#include "AnnotationManager.h"
#include "Brain.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "DisplayGroupAndTabItemInterface.h"
#include "DisplayGroupAndTabItemTreeWidgetItem.h"
#include "DisplayPropertiesAnnotation.h"
#include "DisplayPropertiesSamples.h"
#include "EventGraphicsPaintSoonAllWindows.h"
#include "EventManager.h"
#include "GuiManager.h"

using namespace caret;


    
/**
 * \class caret::DisplayGroupAndTabItemViewController 
 * \brief View controller for display group and tab item hierarchy
 * \ingroup GuiQt
 */

/**
 * Constructor.
 * 
 * @param dataFileType
 *     Type of data file using this view controller.
 * @param browserWindowIndex
 *     The browser window containing this instance.
 * @param parent
 *     Parent of this instance.
 */
DisplayGroupAndTabItemViewController::DisplayGroupAndTabItemViewController(const DataFileTypeEnum::Enum dataFileType,
                                                                           const int32_t browserWindowIndex,
                                                                           QWidget* parent)
: QWidget(parent),
m_dataFileType(dataFileType),
m_browserWindowIndex(browserWindowIndex)
{
    const QString onOffToolTip("<html>"
                               "To select more than one item:<br>"
                               "* For a contiguous selection, click an item "
                               "and then click another item while holding down "
                               "the SHIFT key.  <br>"
                               "* For non-contiguous selection, select items while "
                               "holding down the CTRL key (Command key on Apple)"
                               "</html>");
    m_turnOnSelectedItemsAction = new QAction("On");
    m_turnOnSelectedItemsAction->setToolTip(onOffToolTip);
    m_turnOnSelectedItemsAction->setCheckable(false);
    QObject::connect(m_turnOnSelectedItemsAction, &QAction::triggered,
                     this, &DisplayGroupAndTabItemViewController::turnOnSelectedItemsTriggered);
    QToolButton* turnOnToolButton = new QToolButton();
    turnOnToolButton->setDefaultAction(m_turnOnSelectedItemsAction);
    
    m_turnOffSelectedItemsAction = new QAction("Off");
    m_turnOffSelectedItemsAction->setToolTip(onOffToolTip);
    m_turnOffSelectedItemsAction->setCheckable(false);
    QObject::connect(m_turnOffSelectedItemsAction, &QAction::triggered,
                     this, &DisplayGroupAndTabItemViewController::turnOffSelectedItemsTriggered);
    QToolButton* turnOffToolButton = new QToolButton();
    turnOffToolButton->setDefaultAction(m_turnOffSelectedItemsAction);
    
    m_treeWidget = new QTreeWidget();
    m_treeWidget->setHeaderHidden(true);
    m_treeWidget->setSelectionMode(QTreeWidget::NoSelection);
    
    QObject::connect(m_treeWidget, SIGNAL(itemCollapsed(QTreeWidgetItem*)),
                     this, SLOT(itemWasCollapsed(QTreeWidgetItem*)));
    QObject::connect(m_treeWidget, SIGNAL(itemExpanded(QTreeWidgetItem*)),
                     this, SLOT(itemWasExpanded(QTreeWidgetItem*)));
    QObject::connect(m_treeWidget, SIGNAL(itemChanged(QTreeWidgetItem*, int)),
                     this, SLOT(itemWasChanged(QTreeWidgetItem*, int)));
    QObject::connect(m_treeWidget, SIGNAL(itemSelectionChanged()),
                     this, SLOT(itemsWereSelected()));
    m_treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    QObject::connect(m_treeWidget, &QTreeWidget::customContextMenuRequested,
                     this, &DisplayGroupAndTabItemViewController::displayContextMenu);
    
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->setContentsMargins(0, 0, 0, 0);
    buttonLayout->addWidget(new QLabel("Selected Items: "));
    buttonLayout->addWidget(turnOnToolButton);
    buttonLayout->addSpacing(5);
    buttonLayout->addWidget(turnOffToolButton);
    buttonLayout->addStretch();
    
    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addLayout(buttonLayout);
    layout->addWidget(m_treeWidget, 100);
    
    s_allViewControllers.insert(this);
}

/**
 * Destructor.
 */
DisplayGroupAndTabItemViewController::~DisplayGroupAndTabItemViewController()
{
    s_allViewControllers.erase(this);
}

/**
 * Gets called when items are selected.
 */
void
DisplayGroupAndTabItemViewController::itemsWereSelected()
{
    QList<QTreeWidgetItem*> itemsSelected = m_treeWidget->selectedItems();
    
    
    if ( ! itemsSelected.empty()) {
        
        std::vector<DisplayGroupAndTabItemInterface*> itemInterfacesVector;
        QListIterator<QTreeWidgetItem*> itemsIter(itemsSelected);
        while (itemsIter.hasNext()) {
            QTreeWidgetItem* item = itemsIter.next();
            DisplayGroupAndTabItemTreeWidgetItem* widgetItem = dynamic_cast<DisplayGroupAndTabItemTreeWidgetItem*>(item);
            CaretAssert(widgetItem);
            DisplayGroupAndTabItemInterface* itemInterface = widgetItem->getDisplayGroupAndTabItem();
            CaretAssert(itemInterface);
            if (itemInterface != NULL) {
                itemInterfacesVector.push_back(itemInterface);
            }
        }
        
        if ( ! itemInterfacesVector.empty()) {
            if (m_dataFileType == DataFileTypeEnum::ANNOTATION) {
                processAnnotationDataSelection(itemInterfacesVector);
            }
            else if (m_dataFileType == DataFileTypeEnum::SAMPLES) {
                processAnnotationDataSelection(itemInterfacesVector);
            }
        }
    }
    
    DisplayGroupEnum::Enum displayGroup = DisplayGroupEnum::DISPLAY_GROUP_TAB;
    int32_t tabIndex = -1;
    getDisplayGroupAndTabIndex(displayGroup, tabIndex);
    updateSelectedAndExpandedCheckboxes(displayGroup,
                                        tabIndex);
    
    updateGraphics();
}

/**
 * Display a context sensitive (right-click) menu.
 *
 * @param pos
 *     Position for context menu
 */
void
DisplayGroupAndTabItemViewController::displayContextMenu(const QPoint& pos)
{
    QList<QTreeWidgetItem*> itemsSelected = m_treeWidget->selectedItems();
    
    if (itemsSelected.isEmpty()) {
        return;
    }
    
    QMenu menu(this);
    QAction* onAction = menu.addAction("Turn all selected items ON");
    menu.addAction("Turn all selected items OFF");
    
    QSignalBlocker blocker(m_treeWidget);
    QAction* selectedAction = menu.exec(m_treeWidget->mapToGlobal(pos));
    if (selectedAction == NULL) {
        return;
    }

    const bool newStatus = (selectedAction == onAction);
    setCheckedStatusOfSelectedItems(newStatus);
}

/**
 * Process the selection of annotations.
 *
 * @param interfaceItems
 *     Items that should be annotations !
 */
void
DisplayGroupAndTabItemViewController::processAnnotationDataSelection(const std::vector<DisplayGroupAndTabItemInterface*>& interfaceItems)
{
    std::set<Annotation*> annotationSet;
    
    for (std::vector<DisplayGroupAndTabItemInterface*>::const_iterator iter = interfaceItems.begin();
         iter != interfaceItems.end();
         iter++) {
        Annotation* ann = dynamic_cast<Annotation*>(*iter);
        if (ann != NULL) {
            annotationSet.insert(ann);
        }
        else {
            AnnotationGroup* annGroup = dynamic_cast<AnnotationGroup*>(*iter);
            if (annGroup != NULL) {
                std::vector<Annotation*> groupAnns;
                annGroup->getAllAnnotations(groupAnns);
                
                annotationSet.insert(groupAnns.begin(),
                                     groupAnns.end());
            }
        }
    }
    
    if ( ! annotationSet.empty()) {
        std::vector<Annotation*> selectedAnnotations(annotationSet.begin(),
                                                     annotationSet.end());
        if (m_dataFileType == DataFileTypeEnum::ANNOTATION) {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager(UserInputModeEnum::Enum::ANNOTATIONS);
            annMan->setAnnotationsForEditing(m_browserWindowIndex,
                                             selectedAnnotations);
        }
        else if (m_dataFileType == DataFileTypeEnum::SAMPLES) {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager(UserInputModeEnum::Enum::SAMPLES_EDITING);
            annMan->setAnnotationsForEditing(m_browserWindowIndex,
                                             selectedAnnotations);
        }
    }
}


/**
 * Gets called when an item is collapsed so that its children are not visible.
 *
 * @param item
 *    The QTreeWidgetItem that was collapsed.
 */
void
DisplayGroupAndTabItemViewController::itemWasCollapsed(QTreeWidgetItem* item)
{
    processItemExpanded(item,
                        false);
}

/**
 * Gets called when an item is expaned so that its children are visible.
 *
 * @param item
 *    The QTreeWidgetItem that was expanded.
 */
void
DisplayGroupAndTabItemViewController::itemWasExpanded(QTreeWidgetItem* item)
{
    processItemExpanded(item,
                        true);
}

/**
 * Called when an item is changed (checkbox selected/deselected).
 *
 * @param item
 *    The QTreeWidgetItem that was collapsed.
 * @param column
 *    Ignored.
 */
void
DisplayGroupAndTabItemViewController::itemWasChanged(QTreeWidgetItem* item,
                                                    int /*column*/)
{
    DisplayGroupEnum::Enum displayGroup = DisplayGroupEnum::DISPLAY_GROUP_TAB;
    int32_t tabIndex = -1;
    getDisplayGroupAndTabIndex(displayGroup, tabIndex);
    
    DisplayGroupAndTabItemInterface* dataItem = getDataItem(item);
    
    const Qt::CheckState checkState = item->checkState(DisplayGroupAndTabItemTreeWidgetItem::NAME_COLUMN);
    const TriStateSelectionStatusEnum::Enum itemCheckState = DisplayGroupAndTabItemTreeWidgetItem::fromQCheckState(checkState);
    dataItem->setItemDisplaySelected(displayGroup,
                              tabIndex,
                              itemCheckState);

    updateSelectedAndExpandedCheckboxes(displayGroup,
                                        tabIndex);
    updateSelectedAndExpandedCheckboxesInOtherViewControllers();
    
    updateGraphics();
}

/**
 * Process item expanded or collapsed.
 *
 * @param item
 *     The QTreeWidgetItem that was expanded or collapsed.
 * @param expandedStatus
 *     True if expanded, false if collapsed.
 */
void
DisplayGroupAndTabItemViewController::processItemExpanded(QTreeWidgetItem* item,
                                                          const bool expandedStatus)
{
    DisplayGroupEnum::Enum displayGroup = DisplayGroupEnum::DISPLAY_GROUP_TAB;
    int32_t tabIndex = -1;
    getDisplayGroupAndTabIndex(displayGroup, tabIndex);
    
    DisplayGroupAndTabItemInterface* dataItem = getDataItem(item);
    dataItem->setItemExpanded(displayGroup,
                              tabIndex,
                              expandedStatus);
    updateSelectedAndExpandedCheckboxes(displayGroup,
                                        tabIndex);
    updateSelectedAndExpandedCheckboxesInOtherViewControllers();
    
}

/**
 * Get the data item in the given tree widget item.
 * 
 * @param item
 *      The tree widget item.
 * @return
 *      The data item in the tree widget item.
 */
DisplayGroupAndTabItemInterface*
DisplayGroupAndTabItemViewController::getDataItem(QTreeWidgetItem* item) const
{
    DisplayGroupAndTabItemTreeWidgetItem* treeItem = dynamic_cast<DisplayGroupAndTabItemTreeWidgetItem*>(item);
    CaretAssert(treeItem);
    DisplayGroupAndTabItemInterface* dataItem = treeItem->getDisplayGroupAndTabItem();
    CaretAssert(dataItem);
    return dataItem;
}

/**
 * Get the display group and tab index currently active.
 */
void
DisplayGroupAndTabItemViewController::getDisplayGroupAndTabIndex(DisplayGroupEnum::Enum& displayGroupOut,
                                                                 int32_t& tabIndexOut) const
{
    BrowserTabContent* tabContent = GuiManager::get()->getBrowserTabContentForBrowserWindow(m_browserWindowIndex, false);
    CaretAssert(tabContent);
    tabIndexOut= tabContent->getTabNumber();
    CaretAssert(tabIndexOut >= 0);
    
    if (m_dataFileType == DataFileTypeEnum::ANNOTATION) {
        DisplayPropertiesAnnotation* dpa = GuiManager::get()->getBrain()->getDisplayPropertiesAnnotation();
        displayGroupOut = dpa->getDisplayGroupForTab(tabIndexOut);
    }
    else if (m_dataFileType == DataFileTypeEnum::SAMPLES) {
        DisplayPropertiesSamples* dps(GuiManager::get()->getBrain()->getDisplayPropertiesSamples());
        displayGroupOut = dps->getDisplayGroupForTab(tabIndexOut);
    }
}


/**
 * Update the content.
 *
 * @param contentItemsIn
 *     Items that are displayed.
 * @param displayGroup
 *     The display group.
 * @param tabIndex
 *     Index of the tab
 * @param allowSelectionFlag
 *     Allows selection of items by user clicking items.
 */
void
DisplayGroupAndTabItemViewController::updateContent(std::vector<DisplayGroupAndTabItemInterface*>& contentItemsIn,
                                                    const DisplayGroupEnum::Enum displayGroup,
                                                    const int32_t tabIndex,
                                                    const bool allowSelectionFlag)
{
    if (allowSelectionFlag) {
        m_treeWidget->setSelectionMode(QTreeWidget::ExtendedSelection);
    }
    else {
        m_treeWidget->setSelectionMode(QTreeWidget::NoSelection);
    }
    
    /*
     * Ignore items without children
     */
    std::vector<DisplayGroupAndTabItemInterface*> contentItems;
    for (std::vector<DisplayGroupAndTabItemInterface*>::iterator contIter = contentItemsIn.begin();
         contIter != contentItemsIn.end();
         contIter++) {
        DisplayGroupAndTabItemInterface* item = *contIter;
        if (item->getNumberOfItemChildren() > 0) {
            contentItems.push_back(item);
        }
    }
    
    /*
     * Updating the tree will cause signals so block them until update is done
     */
    m_treeWidget->blockSignals(true);
    
    const int32_t numExistingChildren = m_treeWidget->topLevelItemCount();
    const int32_t numValidChildren    = contentItems.size();
    
    const int32_t numberOfChildrenToAdd = numValidChildren - numExistingChildren;
    for (int32_t i = 0; i < numberOfChildrenToAdd; i++) {
        m_treeWidget->addTopLevelItem(new DisplayGroupAndTabItemTreeWidgetItem(m_browserWindowIndex));
    }
    
    CaretAssert(m_treeWidget->topLevelItemCount() >= numValidChildren);
    
    for (int32_t i = 0; i < numValidChildren; i++) {
        QTreeWidgetItem* treeWidgetChild = m_treeWidget->topLevelItem(i);
        CaretAssert(treeWidgetChild);
        DisplayGroupAndTabItemTreeWidgetItem* dgtChild = dynamic_cast<DisplayGroupAndTabItemTreeWidgetItem*>(treeWidgetChild);
        CaretAssert(dgtChild);
        
            treeWidgetChild->setHidden(false);
            
            CaretAssertVectorIndex(contentItems, i);
            CaretAssert(contentItems[i]);
            DisplayGroupAndTabItemInterface* displayGroupAndTabItem = contentItems[i];
            dgtChild->updateContent(displayGroupAndTabItem,
                                    m_treeWidget,
                                    displayGroup,
                                    tabIndex);
    }
    
    for (int32_t i = (numExistingChildren - 1); i >= numValidChildren; i--) {
        /*
         * Take removes it from the parent but
         * does not destruct it.
         */
        QTreeWidgetItem* item = m_treeWidget->takeTopLevelItem(i);
        delete item;
    }

    updateSelectedAndExpandedCheckboxes(displayGroup,
                                        tabIndex);
    
    /*
     * Allow signals now that updating is done
     */
    m_treeWidget->blockSignals(false);
}

/**
 * Update graphics and, in some circumstances, surface node coloring.
 */
void
DisplayGroupAndTabItemViewController::updateGraphics()
{
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    EventManager::get()->sendEvent(EventGraphicsPaintSoonAllWindows().getPointer());
}

/**
 * Update the selected and expanded checkboxes.
 */
void
DisplayGroupAndTabItemViewController::updateSelectedAndExpandedCheckboxes(const DisplayGroupEnum::Enum displayGroup,
                                                                          const int32_t tabIndex)
{
    m_treeWidget->blockSignals(true);
    
    const int32_t numChildren = m_treeWidget->topLevelItemCount();
    for (int32_t itemIndex = 0; itemIndex < numChildren; itemIndex++) {
        QTreeWidgetItem* treeChild = m_treeWidget->topLevelItem(itemIndex);
        CaretAssert(treeChild);
        
        DisplayGroupAndTabItemTreeWidgetItem* item = dynamic_cast<DisplayGroupAndTabItemTreeWidgetItem*>(treeChild);
        CaretAssert(item);
        
        DisplayGroupAndTabItemInterface* data = item->getDisplayGroupAndTabItem();
        if (data != NULL) {
            item->updateSelectedAndExpandedCheckboxes(displayGroup,
                                                      tabIndex);
        }
    }
    
    const bool itemsSelectedFlag = ( ! m_treeWidget->selectedItems().isEmpty());
    m_turnOnSelectedItemsAction->setEnabled(itemsSelectedFlag);
    m_turnOffSelectedItemsAction->setEnabled(itemsSelectedFlag);
    
    m_treeWidget->blockSignals(false);
}

/**
 * Update the selection and expansion controls in ALL other view controllers.
 * All of them need to be updated since window annotation selection is not
 * affected by the display group and tab selection.
 */
void
DisplayGroupAndTabItemViewController::updateSelectedAndExpandedCheckboxesInOtherViewControllers()
{
    for (std::set<DisplayGroupAndTabItemViewController*>::iterator iter = s_allViewControllers.begin();
         iter != s_allViewControllers.end();
         iter++) {
        DisplayGroupAndTabItemViewController* otherViewController = *iter;
        if (otherViewController != this) {
            if (otherViewController->m_dataFileType == m_dataFileType) {
                DisplayGroupEnum::Enum otherDisplayGroup = DisplayGroupEnum::DISPLAY_GROUP_TAB;
                int32_t otherTabIndex = -1;
                otherViewController->getDisplayGroupAndTabIndex(otherDisplayGroup,
                                                                otherTabIndex);
                otherViewController->updateSelectedAndExpandedCheckboxes(otherDisplayGroup,
                                                                         otherTabIndex);
            }
        }
    }
}

/**
 * Turn on all selected items
 */
void
DisplayGroupAndTabItemViewController::turnOnSelectedItemsTriggered()
{
    setCheckedStatusOfSelectedItems(true);
}

/**
 * Turn off all selected items
 */
void
DisplayGroupAndTabItemViewController::turnOffSelectedItemsTriggered()
{
    setCheckedStatusOfSelectedItems(false);
}

/**
 * Set the checked status of all selected itemsj
 *
 * @param checkedStatus
 *     Checked status
 */
void
DisplayGroupAndTabItemViewController::setCheckedStatusOfSelectedItems(const bool checkedStatus)
{
    QList<QTreeWidgetItem*> itemsSelected = m_treeWidget->selectedItems();
    
    if (itemsSelected.isEmpty()) {
        return;
    }
    
    DisplayGroupEnum::Enum displayGroup = DisplayGroupEnum::DISPLAY_GROUP_TAB;
    int32_t tabIndex = -1;
    getDisplayGroupAndTabIndex(displayGroup, tabIndex);
    
    
    const Qt::CheckState newCheckState = (checkedStatus
                                          ? Qt::Checked
                                          : Qt::Unchecked);
    QListIterator<QTreeWidgetItem*> iter(itemsSelected);
    while (iter.hasNext()) {
        QTreeWidgetItem* item = iter.next();
        DisplayGroupAndTabItemInterface* dataItem = getDataItem(item);
        
        const TriStateSelectionStatusEnum::Enum itemCheckState = DisplayGroupAndTabItemTreeWidgetItem::fromQCheckState(newCheckState);
        dataItem->setItemDisplaySelected(displayGroup,
                                         tabIndex,
                                         itemCheckState);
    }
    
    updateSelectedAndExpandedCheckboxes(displayGroup,
                                        tabIndex);
    updateSelectedAndExpandedCheckboxesInOtherViewControllers();
    
    updateGraphics();}

