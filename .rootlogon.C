{
    TClass::GetClass(
        "ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>");
    TClass::GetClass(
        "ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>");
    // Add palette
    TString path {gSystem->Getenv("ACTROOT")};
    path += "/Core/Data/ColorPalettes/";
    // Get Palette
    path += "managua";
    path += ".txt";
    if(gSystem->AccessPathName(path))
    {
        std::cout << BOLDRED << "Palette named " << name << " does not exist in ActRoot database" << RESET << '\n';
        return;
    }
    // Set!
    gStyle->SetPalette(path);
    gStyle->SetNumberContours(256);
    // Reverse it?
    TColor::InvertPalette();
    // Save palette in files
    TColor::DefinedColors(1);
    // Histogram options
    gStyle->SetOptStat("nmeruoi");
    // Pad margins
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadBottomMargin(0.12);
    // Font type
    gStyle->SetTextFont(22);
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetLegendFont(132);
    gStyle->SetStatFont(132);
    gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetTitleFont(132, "");
    // Title sizes
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetTitleSize(0.06, "");
    // Label size
    gStyle->SetLabelSize(0.055, "XYZ");
    // Ticks
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Legend
    gStyle->SetLegendTextSize(0.05);
    // Axis max digits
    gStyle->SetAxisMaxDigits(4);

    // Selector
    gROOT->ProcessLine("#include \"/media/Data/E796v2/Selector/Selector.h\"");
    gROOT->ProcessLine("#include \"/media/Data/E796v2/Selector/Selector.cxx\"");
}
