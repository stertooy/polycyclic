LoadPackage("polycyclic");
LoadPackage("smallgrp");
TestDirectory(DirectoriesPackageLibrary("polycyclic", "tst"), rec(exitGAP := true));
FORCE_QUIT_GAP(1);
