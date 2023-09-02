all:
	javac -cp . bulkatac2fragments/BamToFragments.java 
	cd .; jar cfm ./bulkatac2fragments.jar Manifest.txt ./bulkatac2fragments/*.class

gitaddall:
	git add */*.java

