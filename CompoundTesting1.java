import chemaxon.formats.MolImporter

String file = 'C:/data/import.sdf'

importer = new MolImporter(file)
importer.grabbingEnabled = true
mol = importer.createMol()

while (importer.read(mol)) {
    String molStr = importer.grabbedMoleculeString
    print "Molecule is $mol"   
}