#include "PCASplit.h"

#include "ActColors.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "TPrincipal.h"

#include <iostream>
#include <map>
#include <vector>

void ActAlgorithm::PCASplit::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
}

void ActAlgorithm::PCASplit::Run()
{
    if(!fIsEnabled)
        return;
    std::cout << "Running PCASplit" << '\n';
    for(auto it = fTPCData->fClusters.begin(); it != fTPCData->fClusters.end(); it++)
    {
        // 1-> Check whether we meet conditions to execute this
        auto& line {it->GetLine()};
        bool isBadFit {line.GetChi2() > 1.5};
        if(!isBadFit)
            continue;

        // 2-> Sort by increasing X
        auto& voxels {it->GetRefToVoxels()};
        std::sort(voxels.begin(), voxels.end(), [](const ActRoot::Voxel& a, const ActRoot::Voxel& b)
                  { return a.GetPosition().X() < b.GetPosition().X(); });

        // 3-> Define window along X dir
        int xWindow {5};
        // Store values
        std::map<int, std::vector<double>> components {};
        for(int i = xWindow; i < voxels.size() - xWindow; i += xWindow)
        {
            TPrincipal pca {3, "ND"};
            double xsum {};
            for(int j = i; j < i + xWindow; j++)
            {
                auto& pos {voxels[j].GetPosition()};
                double x[3] {pos.X(), pos.Y(), pos.Z()};
                pca.AddRow(x);
                xsum += pos.X();
            }

            // Call PCA
            pca.MakePrincipals();
            // Get eigenvalues
            auto& eigenval {*pca.GetEigenValues()};
            components[(int)(xsum / xWindow)] = {eigenval(0), eigenval(1), eigenval(2)};
        }
        if(fIsVerbose)
        {
            std::cout << "PCA for cluster#" << it->GetClusterID() << '\n';
            for(const auto& [x, vals] : components)
            {
                std::cout << " x : " << x << " val0: " << vals[0] << " val1: " << vals[1] << '\n';
            }
        }
    }
}

void ActAlgorithm::PCASplit::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  No parameters yet" << RESET << '\n';
}

extern "C" ActAlgorithm::PCASplit* CreateUserAction()
{
    return new ActAlgorithm::PCASplit;
}
