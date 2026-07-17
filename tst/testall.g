CHECK_IGS@Polycyclic := true;
LoadPackage("polycyclic");
TestDirectory(DirectoriesPackageLibrary("polycyclic", "tst"), rec(exitGAP := true));
FORCE_QUIT_GAP(1);
