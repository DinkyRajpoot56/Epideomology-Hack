import javax.swing.*

def chooser = new JFileChooser()
chooser.setDialogTitle('Select file')
if (chooser.showOpenDialog(null)==JFileChooser.APPROVE_OPTION) {
    fileName = chooser.selectedFile.canonicalPath
    // fileName is then used when creating a MolImporter
} else {
    return 
} 