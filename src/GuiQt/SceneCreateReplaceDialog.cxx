
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

#include <algorithm>

#define __SCENE_CREATE_REPLACE_DIALOG_DECLARE__
#include "SceneCreateReplaceDialog.h"
#undef __SCENE_CREATE_REPLACE_DIALOG_DECLARE__

#include <QAction>
#include <QCheckBox>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QToolButton>
#include <QVBoxLayout>

#include "ApplicationInformation.h"
#include "Brain.h"
#include "BrainBrowserWindow.h"
#include "CaretAssert.h"
#include "CaretPreferences.h"
#include "DataFileException.h"
#include "EventBrowserTabGetAll.h"
#include "EventBrowserTabGetAllViewed.h"
#include "EventImageCapture.h"
#include "EventSceneActive.h"
#include "EventManager.h"
#include "GuiManager.h"
#include "ImageFile.h"
#include "PlainTextStringBuilder.h"
#include "Scene.h"
#include "SceneAttributes.h"
#include "SceneDialog.h"
#include "SceneFile.h"
#include "SceneInfo.h"
#include "SessionManager.h"
#include "WuQMessageBox.h"
#include "WuQtUtilities.h"

using namespace caret;
    
/**
 * \class caret::SceneCreateReplaceDialog 
 * \brief Dialog for creating or replacing a scene.
 * \ingroup GuiQt
 *
 */

#include "Scene.h"

/**
 * Constructor.
 *
 * @param dialogTitle
 *     Title of the dialog.
 * @param parent
 *     Parent on which dialog is displayed.
 * @param sceneFile
 *     Scene file to which scene is added or replaced.
 * @param mode
 *     Scene add/insert/replace mode
 * @param sceneToReplace
 *     If non-NULL, this scene will be replaced.
 */
SceneCreateReplaceDialog::SceneCreateReplaceDialog(const AString& dialogTitle,
                                                   QWidget* parent,
                                                   SceneFile* sceneFile,
                                                   const Mode mode,
                                                   Scene* sceneToInsertOrReplace)
: WuQDialogModal(dialogTitle,
                 parent)
{
    CaretAssert(sceneFile);
    switch (m_mode) {
        case MODE_ADD_NEW_SCENE:
            break;
        case MODE_EDIT_SCENE_INFO:
            CaretAssert(sceneToInsertOrReplace);
            break;
        case MODE_INSERT_NEW_SCENE:
            CaretAssert(sceneToInsertOrReplace);
            break;
        case MODE_REPLACE_SCENE:
            CaretAssert(sceneToInsertOrReplace);
            break;
    }
    
    if ( ! s_previousSelectionsValid) {
        s_previousSelectionsValid = true;
        
        s_previousSelections.m_addAllLoadedFiles          = true;
        s_previousSelections.m_addAllTabs                 = true;
        s_previousSelections.m_addModifiedPaletteSettings = true;
        s_previousSelections.m_addSpecFileNameToScene     = true;
        s_previousSelections.m_cropImage                  = true;
    }
    
    m_sceneFile = sceneFile;
    m_mode      = mode;
    m_sceneToInsertOrReplace = sceneToInsertOrReplace;
    m_sceneThatWasCreated = NULL;
    
    /*
     * Options widget
     */
    QLabel* optionsLabel = new QLabel("Options");
    QWidget* optionsWidget = createSceneOptionsWidget();
    
    /*
     * Create scene information widgets
     */
    QLabel* nameLabel = new QLabel("Name");
    m_nameLineEdit = new QLineEdit();
    
    const QString sceneIdWarningToolTip("The BALSA Scene ID is an identifier created by BALSA.  ALTERING "
                                        "THE BALSA SCENE ID MAY CORRUPT THE SCENE FILE.");
    QLabel* sceneIDLabel = new QLabel("BALSA Scene ID");
    const QString sceneIdLineEditToolTip(sceneIdWarningToolTip
                                         + "  Click the "
                                         "\"Unlock\" button if you must edit the ID.");
    m_balsaSceneIDLineEdit = new QLineEdit();
    m_balsaSceneIDLineEdit->setReadOnly(true);
    m_balsaSceneIDLineEdit->setMaxLength(10);
    WuQtUtilities::setWordWrappedToolTip(m_balsaSceneIDLineEdit,
                                         sceneIdLineEditToolTip);
    
    const QString balsaIdToolTipText("Click this button to unlock (lock) and enable (disable) editing of the "
                                     "BALSA Scene ID.  "
                                     + sceneIdWarningToolTip);
    m_lockUnlockBalseSceneIdToolButton = new QToolButton();
    m_lockUnlockBalseSceneIdToolButton->setText("Unlock");
    WuQtUtilities::setWordWrappedToolTip(m_lockUnlockBalseSceneIdToolButton,
                                         balsaIdToolTipText);
    QObject::connect(m_lockUnlockBalseSceneIdToolButton, &QToolButton::clicked,
                     this, &SceneCreateReplaceDialog::lockUnlockBalseSceneIdToolButtonClicked);
    
    
    QPushButton* addWindowDescriptionPushButton = new QPushButton("Add Window Info");
    QObject::connect(addWindowDescriptionPushButton, &QPushButton::clicked,
                     this, &SceneCreateReplaceDialog::addWindowContentToolButtonClicked);
    
    QLabel* descriptionLabel = new QLabel("Description");
    m_descriptionTextEdit = new QPlainTextEdit();
    
    QHBoxLayout* idLayout = new QHBoxLayout();
    idLayout->setContentsMargins(0, 0, 0, 0);
    idLayout->addWidget(m_balsaSceneIDLineEdit);
    idLayout->addWidget(m_lockUnlockBalseSceneIdToolButton);
    idLayout->addStretch();
    
    const Qt::Alignment labelAlignment = (Qt::AlignLeft | Qt::AlignTop);
    
    int32_t columnCounter = 0;
    const int32_t labelColumn  = columnCounter++;
    const int32_t widgetColumn = columnCounter++;
    QGridLayout* infoGridLayout = new QGridLayout();
    infoGridLayout->setColumnStretch(labelColumn,   0);
    infoGridLayout->setColumnStretch(widgetColumn, 100);
    int32_t rowCounter = 0;
    infoGridLayout->setRowStretch(rowCounter, 0);
    infoGridLayout->addWidget(nameLabel,
                              rowCounter, labelColumn,
                              labelAlignment);
    infoGridLayout->addWidget(m_nameLineEdit,
                              rowCounter, widgetColumn);
    rowCounter++;
    infoGridLayout->addWidget(sceneIDLabel,
                              rowCounter, labelColumn,
                              labelAlignment);
//    infoGridLayout->addWidget(m_balsaSceneIDLineEdit,
//                              rowCounter, widgetColumn);
    infoGridLayout->addLayout(idLayout,
                              rowCounter, widgetColumn);
    rowCounter++;
    infoGridLayout->addWidget(descriptionLabel,
                              rowCounter, labelColumn,
                              labelAlignment);
    infoGridLayout->addWidget(m_descriptionTextEdit,
                              rowCounter, widgetColumn,
                              2, 1);
    rowCounter++;
    infoGridLayout->addWidget(addWindowDescriptionPushButton,
                              rowCounter, labelColumn,
                              labelAlignment);
    infoGridLayout->setRowStretch(rowCounter, 100);
    rowCounter++;
    infoGridLayout->addWidget(optionsLabel,
                              rowCounter, labelColumn,
                              labelAlignment);
    infoGridLayout->addWidget(optionsWidget,
                              rowCounter, widgetColumn);
    rowCounter++;
    
    /*
     * Add the layout to a widget and return the widget.
     */
    QWidget* infoWidget = new QWidget();
    infoWidget->setLayout(infoGridLayout);

    
    QWidget* dialogWidget = new QWidget();
    QVBoxLayout* dialogLayout = new QVBoxLayout(dialogWidget);
    dialogLayout->addWidget(infoWidget,
                            100);
    
    setCentralWidget(dialogWidget,
                     WuQDialog::SCROLL_AREA_NEVER);

    PlainTextStringBuilder windowDescriptionBuilder;
    std::vector<BrainBrowserWindow*> windows = GuiManager::get()->getAllOpenBrainBrowserWindows();
    for (std::vector<BrainBrowserWindow*>::iterator iter = windows.begin();
         iter != windows.end();
         iter++) {
        BrainBrowserWindow* window = *iter;
        window->getDescriptionOfContent(windowDescriptionBuilder);
        windowDescriptionBuilder.addLine("");
    }
    const AString windowDescription = windowDescriptionBuilder.getText().trimmed();
    
    QDateTime dateTime = QDateTime::currentDateTime();
    const QString dateTimeText = dateTime.toString("dd MMM yyyy hh:mm:ss");
    const QString commitName   = (ApplicationInformation().getCommit());
    const QString dataTimeCommitText(dateTimeText + "\n" + commitName);
    
    const AString defaultNewSceneName = ("New Scene " + AString::number(sceneFile->getNumberOfScenes() + 1));
    m_nameLineEdit->setText(defaultNewSceneName);
    
    switch (m_mode) {
        case MODE_ADD_NEW_SCENE:
        case MODE_INSERT_NEW_SCENE:
            m_sceneWindowDescription = ("Created on " + dataTimeCommitText + "\n");
            m_sceneWindowDescription.appendWithNewLine(windowDescription);
            break;
        case MODE_EDIT_SCENE_INFO:
            m_nameLineEdit->setText(sceneToInsertOrReplace->getName());
            m_descriptionTextEdit->setPlainText(sceneToInsertOrReplace->getDescription());
            m_balsaSceneIDLineEdit->setText(sceneToInsertOrReplace->getBalsaSceneID());
            addWindowDescriptionPushButton->setHidden(true);
            optionsLabel->setHidden(true);
            optionsWidget->setHidden(true);
            break;
        case MODE_REPLACE_SCENE:
            m_nameLineEdit->setText(sceneToInsertOrReplace->getName());
            m_balsaSceneIDLineEdit->setText(sceneToInsertOrReplace->getBalsaSceneID());
            m_sceneWindowDescription = ("Replaced on " + dataTimeCommitText + "\n");
            m_sceneWindowDescription.appendWithNewLine(windowDescription);
            m_sceneWindowDescription.appendWithNewLine(" ");
            m_sceneWindowDescription.appendWithNewLine(sceneToInsertOrReplace->getDescription());
            m_descriptionTextEdit->setPlainText(sceneToInsertOrReplace->getDescription());
            break;
    }
    
    
    m_nameLineEdit->setFocus();
    m_nameLineEdit->selectAll();
    
    setMinimumWidth(600);
    setMinimumHeight(300);
    
    setSaveWindowPositionForNextTime("SceneCreateDialog");
}

/**
 * Destructor.
 */
SceneCreateReplaceDialog::~SceneCreateReplaceDialog()
{
}

/**
 * Static method that creates a dialog for creating a new scene and 
 * returns the scene that was created or NULL if scene was not created.
 *
 * @param parent
 *     Parent widget on which dialog is displayed.
 * @param sceneFile
 *     Scene file to which new scene is added.
 * @return
 *     Scene that was created or NULL if user cancelled or there was an error.
 */
Scene*
SceneCreateReplaceDialog::createNewScene(QWidget* parent,
                                         SceneFile* sceneFile)
{
    SceneCreateReplaceDialog dialog("Create New Scene",
                                    parent,
                                    sceneFile,
                                    MODE_ADD_NEW_SCENE,
                                    NULL);
    dialog.exec();
    
    Scene* scene = dialog.m_sceneThatWasCreated;
    return scene;
}

/**
 * Called when the scene ID lock unlock button is clicked
 */
void
SceneCreateReplaceDialog::lockUnlockBalseSceneIdToolButtonClicked()
{
    if (m_balsaSceneIDLineEdit->isReadOnly()) {
        m_balsaSceneIDLineEdit->setReadOnly(false);
        const AString text("WARNING: Are you sure you want to edit the BALSA Scene ID?");
        const AString infoText("Changing the BALSA Scene ID may corrupt the Scene File.  Editing of the "
                               "Scene ID should only be performed by those with expert knowledge of the BALSA "
                               "Database.");
        if (WuQMessageBox::warningYesNo(m_lockUnlockBalseSceneIdToolButton,
                                        text,
                                        infoText,
                                        WuQMessageBox::DefaultButtonYesNo::NO)) {
            m_lockUnlockBalseSceneIdToolButton->setText("Lock");
        }
    }
    else {
        m_balsaSceneIDLineEdit->setReadOnly(true);
        m_lockUnlockBalseSceneIdToolButton->setText("Unlock");
    }
}

/**
 * Static method that creates a dialog for editing a scene's info
 *
 * @param parent
 *     Parent widget on which dialog is displayed.
 * @param sceneFile
 *     Scene file to which new scene is added.
 * @param scene
 *     Scene that is edited
 * @return
 *     Scene that was created or NULL if user cancelled or there was an error.
 */
void
SceneCreateReplaceDialog::editSceneInfo(QWidget* parent,
                                        SceneFile* sceneFile,
                                        Scene* scene)
{
    SceneCreateReplaceDialog dialog("Edit Scene",
                                    parent,
                                    sceneFile,
                                    MODE_EDIT_SCENE_INFO,
                                    scene);
    dialog.exec();
}

/**
 * Static method that creates a dialog for creating a new scene and
 * returns the scene that was created or NULL if scene was not created.
 *
 * @param parent
 *     Parent widget on which dialog is displayed.
 * @param sceneFile
 *     Scene file to which new scene is added.
 * @param insertBeforeScene
 *     Insert the newly created scene BEFORE this scene.
 * @return
 *     Scene that was created or NULL if user cancelled or there was an error.
 */
Scene*
SceneCreateReplaceDialog::createNewSceneInsertBeforeScene(QWidget* parent,
                                             SceneFile* sceneFile,
                                             const Scene* insertBeforeScene)
{
    SceneCreateReplaceDialog dialog("Insert New Scene",
                                    parent,
                                    sceneFile,
                                    MODE_INSERT_NEW_SCENE,
                                    const_cast<Scene*>(insertBeforeScene));
    dialog.exec();
    
    Scene* scene = dialog.m_sceneThatWasCreated;
    return scene;
}


/**
 * Static method that creates a dialog for replacing and existing scene and
 * returns the scene that was created or NULL if scene was not created.
 *
 * @param parent
 *     Parent widget on which dialog is displayed.
 * @param sceneFile
 *     File in which the given scene exists and will be replaced.
 * @param sceneToReplace
 *     Scene that is being replaced.  If the user presses the OK button
 *     to replace this scene it will be destroyed so the pointer must not
 *     be deferenced at any time after calling this method.
 * @return
 *     Scene that was created or NULL if user cancelled or there was an error.
 */
Scene*
SceneCreateReplaceDialog::replaceExistingScene(QWidget* parent,
                                               SceneFile* sceneFile,
                                               Scene* sceneToReplace)
{
    const AString title = ("Replace Scene: "
                           + sceneToReplace->getName());
    
    SceneCreateReplaceDialog dialog(title,
                                    parent,
                                    sceneFile,
                                    MODE_REPLACE_SCENE,
                                    sceneToReplace);
    
    /*
     * Run the dialog.
     * If user cancels, "dialog.m_sceneThatWasCreated" will be NULL.
     */
    dialog.exec();
    Scene* scene = dialog.m_sceneThatWasCreated;
    
    return scene;
}

/**
 * Add an image to the scene and also info about Workbench.
 *
 * @param scene
 *    Scene to which image is added.
 * @param cropImageFlag
 *   If true, crop the image
 * @param errorMessageOut
 *   Output with any error messages
 */
void
SceneCreateReplaceDialog::addImageAndWorkbenchInfoToScene(Scene* scene,
                                                          const bool cropImageFlag,
                                                          AString& errorMessageOut)
{
    errorMessageOut.clear();
    
    CaretAssert(scene);

    uint8_t backgroundColor[3] = { 0, 0, 0 };
    bool backgroundColorValid = false;

    /*
     * Capture an image of each window
     */
    std::vector<ImageFile*> imageFiles;
    std::vector<BrainBrowserWindow*> windows = GuiManager::get()->getAllOpenBrainBrowserWindows();
    for (std::vector<BrainBrowserWindow*>::iterator iter = windows.begin();
         iter != windows.end();
         iter++) {
        BrainBrowserWindow* bbw = *iter;
        const int32_t browserWindowIndex = bbw->getBrowserWindowIndex();
        
        EventImageCapture imageCaptureEvent(browserWindowIndex);
        EventManager::get()->sendEvent(imageCaptureEvent.getPointer());
        
        if (imageCaptureEvent.getEventProcessCount() > 0) {
            if (imageCaptureEvent.isError()) {
                errorMessageOut.appendWithNewLine(imageCaptureEvent.getErrorMessage());
            }
            else {
                imageFiles.push_back(new ImageFile(imageCaptureEvent.getCapturedImage()));
                if ( ! backgroundColorValid) {
                    imageCaptureEvent.getBackgroundColor(backgroundColor);
                    backgroundColorValid = true;
                }
            }
        }
    }
    
    /*
     * Assemble images of each window into a single image
     * and add it to the scene.  Use one image per row
     * since the images are limited in horizontal space
     * when shown in the listing of scenes.
     */
    if ( ! imageFiles.empty()) {
        try {
            const int32_t numImagesPerRow = 1;
            ImageFile compositeImageFile;
            compositeImageFile.combinePreservingAspectAndFillIfNeeded(imageFiles,
                                                                      numImagesPerRow,
                                                                      backgroundColor);
            
            if (backgroundColorValid) {
                if (cropImageFlag) {
                    const int marginSize = 5;
                    compositeImageFile.cropImageRemoveBackground(marginSize,
                                                                 backgroundColor);
                }
            }
            
            const int MAXIMUM_IMAGE_WIDTH = 1024;
            compositeImageFile.resizeToMaximumWidth(MAXIMUM_IMAGE_WIDTH);
            
            const AString PREFERRED_IMAGE_FORMAT = "png";
            
            QByteArray byteArray;
            compositeImageFile.getImageInByteArray(byteArray,
                                                   PREFERRED_IMAGE_FORMAT);
            
            scene->getSceneInfo()->setImageBytes(byteArray,
                                                 PREFERRED_IMAGE_FORMAT);
        }
        catch (const DataFileException& dfe) {
            errorMessageOut.appendWithNewLine((dfe.whatString()
                                               + "\n\nEven though image failed, scene was created."));
        }
    }
    
    /*
     * Free memory from the image files.
     */
    for (std::vector<ImageFile*>::iterator iter = imageFiles.begin();
         iter != imageFiles.end();
         iter++) {
        delete *iter;
    }
    
    scene->getSceneInfo()->addWorkbenchVersionInfoToSceneMetaData();
}

/**
 * Create an image for the loaded scene.
 *
 * @param imageOut
 *     Output image of the scene.
 * @param cropImageFlag
 *     If true, crop the image
 * @param errorMessageOut
 *     Contains error information if image was not created.
 * @return
 *     True if output image is valid, else false.
 */
bool
SceneCreateReplaceDialog::createSceneImage(QImage& imageOut,
                                           const bool cropImageFlag,
                                           AString& errorMessageOut)
{
    bool validImageFlag = false;
    imageOut = QImage();
    errorMessageOut.clear();
    
    uint8_t backgroundColor[3] = { 0, 0, 0 };
    bool backgroundColorValid = false;
    
    /*
     * Capture an image of each window
     */
    std::vector<ImageFile*> imageFiles;
    std::vector<BrainBrowserWindow*> windows = GuiManager::get()->getAllOpenBrainBrowserWindows();
    for (std::vector<BrainBrowserWindow*>::iterator iter = windows.begin();
         iter != windows.end();
         iter++) {
        BrainBrowserWindow* bbw = *iter;
        const int32_t browserWindowIndex = bbw->getBrowserWindowIndex();
        
        EventImageCapture imageCaptureEvent(browserWindowIndex);
        EventManager::get()->sendEvent(imageCaptureEvent.getPointer());
        
        if (imageCaptureEvent.getEventProcessCount() > 0) {
            if (imageCaptureEvent.isError()) {
                errorMessageOut.appendWithNewLine(imageCaptureEvent.getErrorMessage());
            }
            else {
                imageFiles.push_back(new ImageFile(imageCaptureEvent.getCapturedImage()));
                if ( ! backgroundColorValid) {
                    imageCaptureEvent.getBackgroundColor(backgroundColor);
                    backgroundColorValid = true;
                }
            }
        }
    }
    
    /*
     * Assemble images of each window into a single image
     * and add it to the scene.  Use one image per row
     * since the images are limited in horizontal space
     * when shown in the listing of scenes.
     */
    if ( ! imageFiles.empty()) {
        try {
            const int32_t numImagesPerRow = 1;
            ImageFile compositeImageFile;
            compositeImageFile.combinePreservingAspectAndFillIfNeeded(imageFiles,
                                                                      numImagesPerRow,
                                                                      backgroundColor);
            
            if (backgroundColorValid) {
                if (cropImageFlag) {
                    const int marginSize = 5;
                    compositeImageFile.cropImageRemoveBackground(marginSize,
                                                                 backgroundColor);
                }
            }
            
            const int MAXIMUM_IMAGE_WIDTH = 1024;
            compositeImageFile.resizeToMaximumWidth(MAXIMUM_IMAGE_WIDTH);
            
            const AString PREFERRED_IMAGE_FORMAT = "png";
            
            QByteArray byteArray;
            compositeImageFile.getImageInByteArray(byteArray,
                                                   PREFERRED_IMAGE_FORMAT);
            
            imageOut = *compositeImageFile.getAsQImage();
            validImageFlag = true;
        }
        catch (const DataFileException&) {
            errorMessageOut.appendWithNewLine("Even though image generation failed, scene was created.");
        }
    }
    
    /*
     * Free memory from the image files.
     */
    for (std::vector<ImageFile*>::iterator iter = imageFiles.begin();
         iter != imageFiles.end();
         iter++) {
        delete *iter;
    }
    
    return validImageFlag;
}

/**
 * @return Widget containing the scene options widgets.
 */
QWidget*
SceneCreateReplaceDialog::createSceneOptionsWidget()
{
    /*
     * Create scene options widgets
     */
    m_addSpecFileNameToSceneCheckBox = new QCheckBox("Add name of spec file to scene");
    m_addSpecFileNameToSceneCheckBox->setChecked(s_previousSelections.m_addSpecFileNameToScene);
    WuQtUtilities::setWordWrappedToolTip(m_addSpecFileNameToSceneCheckBox,
                                         "Include name of spec file in the scene");
    
    m_addAllTabsCheckBox = new QCheckBox("Add all tabs to scene");
    m_addAllTabsCheckBox->setChecked(s_previousSelections.m_addAllTabs);
    WuQtUtilities::setWordWrappedToolTip(m_addAllTabsCheckBox,
                                         "Add all tabs to the scene.  When this option is selected, "
                                         "the scene will be larger and require additional time to "
                                         "load.  If NOT selected, only the selected tab in each "
                                         "window is saved to the scene.");
    
    m_addAllLoadedFilesCheckBox = new QCheckBox("Add all loaded files to scene");
    m_addAllLoadedFilesCheckBox->setChecked(s_previousSelections.m_addAllLoadedFiles);
    WuQtUtilities::setWordWrappedToolTip(m_addAllLoadedFilesCheckBox,
                                         "Add all loaded files to scene.  When this option is selected, "
                                         "the scene may require additional time to load as file that "
                                         "play no role in reproducing the scene will be loaded.  If NOT "
                                         "selected, the scene may load more quickly.");
    
    m_addModifiedPaletteSettingsCheckBox = new QCheckBox("Add modified palette color mapping to scene");
    m_addModifiedPaletteSettingsCheckBox->setChecked(s_previousSelections.m_addModifiedPaletteSettings);
    WuQtUtilities::setWordWrappedToolTip(m_addModifiedPaletteSettingsCheckBox,
                                         "The palette color mapping is saved within each data files that maps "
                                         "its data to brainordinates.  However, there are instances in which "
                                         "the user wants the scene to display the data with palette color mapping "
                                         "that is different from that in the file.  If this option is "
                                         "selected, modified palettes color mapping will be saved to the scene "
                                         "and the data files with modified palette color mapping do not need "
                                         "to be saved.");
    
    m_cropImageCheckBox = new QCheckBox("Crop Image");
    m_cropImageCheckBox->setChecked(s_previousSelections.m_cropImage);
    WuQtUtilities::setWordWrappedToolTip(m_cropImageCheckBox,
                                         "Crop the image by removing background pixels on the sides of the "
                                         "brain models.  "
                                         "Disabling of image cropping may be useful when "
                                         "debugging scenes where the content occupies a "
                                         "portion of the viewing region (data may be "
                                         "panned/zoomed");
    
    /*
     * Layout for scene options widgets
     */
    QVBoxLayout* optionsLayout = new QVBoxLayout();
    optionsLayout->addWidget(m_addSpecFileNameToSceneCheckBox);
    optionsLayout->addWidget(m_addAllTabsCheckBox);
    optionsLayout->addWidget(m_addAllLoadedFilesCheckBox);
    optionsLayout->addWidget(m_addModifiedPaletteSettingsCheckBox);
    optionsLayout->addWidget(m_cropImageCheckBox);
    
    /*
     * Add the layout to a widget and return the widget.
     */
    QFrame* optionsWidget = new QFrame();
    optionsWidget->setFrameStyle(QFrame::Box
                         | QFrame::Plain);
    optionsWidget->setLineWidth(1);
//    QWidget* optionsWidget = new QWidget();
    optionsWidget->setLayout(optionsLayout);
    return optionsWidget;
}

/**
 * Gets called if the user presses the OK button.
 */
void
SceneCreateReplaceDialog::okButtonClicked()
{
    s_previousSelections.m_addAllLoadedFiles          = m_addAllLoadedFilesCheckBox->isChecked();
    s_previousSelections.m_addAllTabs                 = m_addAllTabsCheckBox->isChecked();
    s_previousSelections.m_addModifiedPaletteSettings = m_addModifiedPaletteSettingsCheckBox->isChecked();
    s_previousSelections.m_addSpecFileNameToScene     = m_addSpecFileNameToSceneCheckBox->isChecked();
    s_previousSelections.m_cropImage                  = m_cropImageCheckBox->isChecked();
    
    const AString newSceneName = m_nameLineEdit->text();
    
    AString errorMessage;
    if (newSceneName.isEmpty()) {
        errorMessage = "Scene Name is empty.";
    }
    else {
        const Scene* sceneWithName = m_sceneFile->getSceneWithName(newSceneName);
        if (sceneWithName != NULL) {
            bool nameErrorFlag = true;
            switch (m_mode) {
                case MODE_ADD_NEW_SCENE:
                    break;
                case MODE_EDIT_SCENE_INFO:
                    if (m_sceneToInsertOrReplace == sceneWithName) {
                        nameErrorFlag = false;
                    }
                    break;
                case MODE_INSERT_NEW_SCENE:
                    break;
                case MODE_REPLACE_SCENE:
                    if (m_sceneToInsertOrReplace == sceneWithName) {
                        nameErrorFlag = false;
                    }
                    break;
            }
            
            if (nameErrorFlag) {
                errorMessage = ("An existing scene uses the name \""
                                + newSceneName
                                + "\".  Scene names must be unique.");
            }
        }
    }
    
    if ( ! errorMessage.isEmpty()) {
        WuQMessageBox::errorOk(this,
                               errorMessage);
        return;
    }

    switch (m_mode) {
        case MODE_ADD_NEW_SCENE:
            break;
        case MODE_EDIT_SCENE_INFO:
        {
            m_sceneToInsertOrReplace->setName(newSceneName);
            m_sceneToInsertOrReplace->setDescription(m_descriptionTextEdit->toPlainText());
            m_sceneToInsertOrReplace->setBalsaSceneID(m_balsaSceneIDLineEdit->text().trimmed());
            WuQDialogModal::okButtonClicked();
            return;
        }
            break;
        case MODE_INSERT_NEW_SCENE:
            break;
        case MODE_REPLACE_SCENE:
            break;
    }

    if ( ! s_previousSelections.m_addModifiedPaletteSettings) {
        if ( ! SceneDialog::checkForModifiedFiles(GuiManager::TEST_FOR_MODIFIED_FILES_PALETTE_ONLY_MODE_FOR_SCENE_ADD,
                                                  this)) {
            /*
             * Add modified palettes to scene is off but
             * there are modified palettes and user has
             * chose to not create the scene
             */
            return;
        }
    }
    
    Scene* newScene = new Scene(SceneTypeEnum::SCENE_TYPE_FULL);
    Scene::setSceneBeingCreated(newScene);
    newScene->setName(newSceneName);
    newScene->setDescription(m_descriptionTextEdit->toPlainText());
    newScene->setBalsaSceneID(m_balsaSceneIDLineEdit->text().trimmed());

    const std::vector<int32_t> windowIndices = GuiManager::get()->getAllOpenBrainBrowserWindowIndices();
    
    /*
     * Get all browser tabs and only save transformations for tabs
     * that are valid.
     */
    std::vector<int32_t> tabIndices;
    if (m_addAllTabsCheckBox->isChecked()) {
        EventBrowserTabGetAll getAllTabs;
        EventManager::get()->sendEvent(getAllTabs.getPointer());
        tabIndices = getAllTabs.getBrowserTabIndices();
    }
    else {
        EventBrowserTabGetAllViewed getViewedTabs;
        EventManager::get()->sendEvent(getViewedTabs.getPointer());
        tabIndices = getViewedTabs.getViewdedBrowserTabIndices();
    }
    std::sort(tabIndices.begin(),
              tabIndices.end());
    
    SceneAttributes* sceneAttributes = newScene->getAttributes();
    sceneAttributes->setSceneFileName(m_sceneFile->getFileName());
    sceneAttributes->setSceneName(newSceneName);
    sceneAttributes->setIndicesOfTabsAndWindowsForSavingToScene(tabIndices,
                                                                windowIndices);
    sceneAttributes->setSpecFileNameSavedToScene(m_addSpecFileNameToSceneCheckBox->isChecked());
    sceneAttributes->setAllLoadedFilesSavedToScene(m_addAllLoadedFilesCheckBox->isChecked());
    sceneAttributes->setModifiedPaletteSettingsSavedToScene(m_addModifiedPaletteSettingsCheckBox->isChecked());
    
    newScene->addClass(GuiManager::get()->saveToScene(sceneAttributes,
                                                      "guiManager"));
    
    AString imageErrorMessage;
    addImageAndWorkbenchInfoToScene(newScene,
                                    s_previousSelections.m_cropImage,
                                    imageErrorMessage);
    if ( ! imageErrorMessage.isEmpty()) {
        WuQMessageBox::errorOk(this,
                               imageErrorMessage);
    }
    
    /*
     * Copy macros from active scene to the new scene
     */
    EventSceneActive activeSceneEvent(EventSceneActive::MODE_GET);
    EventManager::get()->sendEvent(activeSceneEvent.getPointer());
    const Scene* activeScene = activeSceneEvent.getScene();
    if (activeScene != NULL) {
        newScene->copyMacrosFromScene(activeScene);
    }
    
    switch (m_mode) {
        case MODE_ADD_NEW_SCENE:
            m_sceneFile->addScene(newScene);
            break;
        case MODE_EDIT_SCENE_INFO:
            CaretAssert(0);
            break;
        case MODE_INSERT_NEW_SCENE:
            m_sceneFile->insertScene(newScene,
                                     m_sceneToInsertOrReplace);
            break;
        case MODE_REPLACE_SCENE:
            m_sceneFile->replaceScene(newScene,
                                      m_sceneToInsertOrReplace);
            break;
    }
    
    m_sceneThatWasCreated = newScene;
    
    Scene::setSceneBeingCreated(NULL);
    
    WuQDialogModal::okButtonClicked();
}

/**
 * Called when add window content button clicked
 */
void
SceneCreateReplaceDialog::addWindowContentToolButtonClicked()
{
    AString txt = m_descriptionTextEdit->document()->toPlainText();
    txt.appendWithNewLine(m_sceneWindowDescription);
    m_descriptionTextEdit->setPlainText(txt);
}


