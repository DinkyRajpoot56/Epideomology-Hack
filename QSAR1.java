/**
 * Returns the type of the descriptor as defined in the descriptor dictionary.
 * 
 * The method will look for the identifier specified by the user in the QSAR descriptor
 * dictionary. If a corresponding entry is found, first child element that is called
 * "isClassifiedAs" is returned. Note that the OWL descriptor spec allows both the class of
 * descriptor (electronic, topological etc) as well as the type of descriptor (molecular, atomic)
 * to be specified in an "isClassifiedAs" element. Thus we ignore any such element that
 * indicates the descriptors class.
 * 
 * The method assumes that any descriptor entry will have only one "isClassifiedAs" entry describing
 * the descriptors type.
 * 
 * The descriptor can be identified it DescriptorSpecification object
 *
 * @param descriptorSpecification A DescriptorSpecification object
 * @return he type of the descriptor as stored in the dictionary, null if no entry is found matching
 *         the supplied identifier
 */
public String getDictionaryType(IImplementationSpecification descriptorSpecification) {
  return getDictionaryType(descriptorSpecification.getSpecificationReference());
}
origin: cdk/cdk
DescriptorEngineTest.testDictionaryType()
@Test
public void testDictionaryType() {
  DescriptorEngine engine = new DescriptorEngine(IMolecularDescriptor.class,
      DefaultChemObjectBuilder.getInstance());
  String className = "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor";
  DescriptorSpecification specRef = new DescriptorSpecification(
      "http://www.blueobelisk.org/ontologies/chemoinformatics-algorithms/#zagrebIndex", this.getClass()
          .getName(), "The Chemistry Development Kit");
  Assert.assertEquals("molecularDescriptor", engine.getDictionaryType(className));
  Assert.assertEquals("molecularDescriptor", engine.getDictionaryType(specRef));
}
getDictionaryType() of org.openscience.cdk.qsar.DescriptorEngine
Returns the type of the descriptor as defined in the descriptor dictionary. The method will look for the identifier specified by the user in the QSAR descriptor dictionary. If a corresponding entry is found, first child element that is called "isClassifiedAs" is returned. Note that the OWL descriptor spec allows both the class of descriptor (electronic, topological etc) as well as the type of descriptor (molecular, atomic) to be specified in an "isClassifiedAs" element. Thus we ignore any such element that indicates the descriptors class. The method assumes that any descriptor entry will have only one "isClassifiedAs" entry describing the descriptors type. The descriptor can be identified it DescriptorSpecification object

org.openscience.cdk.qsar
DescriptorEngine
getDictionaryType
Javadoc
Show more
Returns the type of the descriptor as defined in the descriptor dictionary. The method will look for the identifier specified by the user in the QSAR descriptor dictionary. If a corresponding entry is found, first child element that is called "isClassifiedAs" is returned. Note that the OWL descriptor spec allows both the class of descriptor (electronic, topological etc) as well as the type of descriptor (molecular, atomic) to be specified in an "isClassifiedAs" element. Thus we ignore any such element that indicates the descriptors class. The method assumes that any descriptor entry will have only one "isClassifiedAs" entry describing the descriptors type. The descriptor can be identified either by the name of the class implementing the descriptor or else the specification reference value of the descriptor which can be obtained from an instance of the descriptor class.
Popular methods of DescriptorEngine
getDictionaryClass
Returns the class(es) of the descriptor as defined in the descriptor dictionary. The method will loo
<init>
Instantiates the DescriptorEngine. This constructor instantiates the engine but does not perform any
getAvailableDictionaryClasses
Get the all the unique dictionary classes that the descriptors belong to.
getDescriptorClassNames
Returns a list containing the names of the classes implementing the descriptors.
getDescriptorInstances
Returns a List containing the instantiated descriptor classes.
getDescriptorSpecifications
Returns the DescriptorSpecification objects for all available descriptors.
getDictionaryDefinition
Gets the definition of the descriptor. All descriptors in the descriptor dictioanry will have a defi
getDictionaryTitle
Gets the label (title) of the descriptor.
getSpecRef
initializeSpecifications
instantiate
instantiateDescriptors
instantiate,instantiateDescriptors,process
Popular in Java
Reactive rest calls using spring rest template
getContentResolver (Context)
getSupportFragmentManager (FragmentActivity)
scheduleAtFixedRate (Timer)
Menu (java.awt)
Proxy (java.net)
This class represents proxy server settings. A created instance of Proxy stores a type and an addres
Path (java.nio.file)
Collections (java.util)
This class consists exclusively of static methods that operate on or return collections. It contains
PriorityQueue (java.util)
A PriorityQueue holds elements on a priority heap, which orders the elements according to their natu
Annotation (javassist.bytecode.annotation)
The annotation structure.An instance of this class is returned bygetAnnotations() in AnnotationsAttr
From CI to AI: The AI layer in your organization
Get Tabnine AI autocomplete in your IDE
