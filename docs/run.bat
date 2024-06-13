@REM echo y| del -r source\apidocs\*

@REM sphinx-apidoc -f -o source/apidocs ../gamut/test

sphinx-apidoc -f -o source/apidoc ../gamut

rm -rf build
sphinx-autobuild source build
