LoadPackage("polycyclic");
CHECK_CENT@      := true;
CHECK_IGS@       := true;
CHECK_INTNORM@   := true;
CHECK_INTSTAB@   := true;
CHECK_NORM@      := true;
CHECK_SCHUR_PCP@ := true;
TestDirectory(DirectoriesPackageLibrary("polycyclic", "tst"), rec(exitGAP := true));
FORCE_QUIT_GAP(1);
