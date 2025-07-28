
//=========================================================================================================================================
// Este es el codigo que hace el algoritmo mio para todos los eventos 
//=========================================================================================================================================


#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <TCanvas.h>
#include <TH1F.h>
#include "TFile.h"
#include "ROOT/RVec.hxx"

// ============================
// Multi-threading
//ROOT can benefit from multi-threading by being able to do work in parallel. If parallelism is managed internally by ROOT, we talk about 
//implicit multi-threading. On the other hand, if the application itself is in charge of creating threads and distributing work among them, 
//the multi-threading is explicit. Below we describe how to use ROOT in either case.
// ============================
ROOT::EnableImplicitMT();  // Activar procesamiento paralelo

// En este caso lo que hace el multireading es decirle al computador, ey haz todo este codigo pero repartido entre cada nucleo del CPU, 
// esto lo hace implicito, por ende no tengo que indicar como operar en cada nucleo, el lo hace completamente solo


// Funcion para calcular distancia en el espacio eta-phi con la funcion de DeltaPhi del TlorentzVector

// ====================================================
// Aqui usa el calculo de ΔR y encuentra los jets mas cercanos para jet_0 y jet_1 sin que se calcule el ΔR entre ellos 
// ====================================================

std::pair<int, int> findClosestJets(const ROOT::VecOps::RVec<float>& phi,
                                    const ROOT::VecOps::RVec<float>& eta,
                                    const ROOT::VecOps::RVec<float>& pt,
                                    const ROOT::VecOps::RVec<float>& mass) {
    
    int idx0 = -1;
    double minDeltaR0 = 999.0;
    TLorentzVector jet0, jet_i;
    jet0.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
    
    for (size_t i = 2; i < phi.size(); ++i) {
        jet_i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        double dr = jet0.DeltaR(jet_i);  // Cambio directo aquí
        if (dr < minDeltaR0) {
            minDeltaR0 = dr;
            idx0 = i;
        }
    }

    int idx1 = -1;
    double minDeltaR1 = 999.0;
    TLorentzVector jet1;
    jet1.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
    
    for (size_t i = 2; i < phi.size(); ++i) {
        if (i == static_cast<size_t>(idx0)) continue;
        jet_i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        double dr = jet1.DeltaR(jet_i);  // Cambio directo aquí
        if (dr < minDeltaR1) {
            minDeltaR1 = dr;
            idx1 = i;
        }
    }
    
    if (idx0 == idx1 && idx0 != -1) {
        jet_i.SetPtEtaPhiM(pt[idx0], eta[idx0], phi[idx0], mass[idx0]);
        double dr_leader = jet0.DeltaR(jet_i);  // Cambio directo aquí
        double dr_subl = jet1.DeltaR(jet_i);    // Cambio directo aquí
        
        if (dr_leader < dr_subl) {
            // Jet x se asigna al líder. Sublíder busca otro jet.
            minDeltaR1 = 999.0;
            for (size_t i = 2; i < phi.size(); ++i) {
                if (i == static_cast<size_t>(idx0)) continue;
                jet_i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
                double dr = jet1.DeltaR(jet_i);  // Cambio directo aquí
                if (dr < minDeltaR1) {
                    minDeltaR1 = dr;
                    idx1 = i;
                }
            }
        } else {
            // Jet x se asigna al sublíder. Líder busca otro jet.
            minDeltaR0 = 999.0;
            for (size_t i = 2; i < phi.size(); ++i) {
                if (i == static_cast<size_t>(idx1)) continue;
                jet_i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
                double dr = jet0.DeltaR(jet_i);  // Cambio directo aquí
                if (dr < minDeltaR0) {
                    minDeltaR0 = dr;
                    idx0 = i;
                }
            }
        }
    }

    return {idx0, idx1};
}
// ===================================================
// Orden descendente en pT de los Jets, este orden lo hace para todas las variables acorde a como se haya ordenado por pT
// ===================================================

// Aqui no se definen como const porque se va a modificar en este caso la posicion que ocupanban anteriormente, acorde pT pues
void sortJetsByPt(ROOT::VecOps::RVec<float>& pt,
                  ROOT::VecOps::RVec<float>& eta,
                  ROOT::VecOps::RVec<float>& phi,
                  ROOT::VecOps::RVec<float>& mass) {
                  
                  
    // Se crea un vector con los indices originales, es decir con las posiciones originales pues
    std::vector<size_t> indices(pt.size());
     // esta es una funcion de C++. este es un vector que se va llenando poco acorde a como yo le indice y a cuantas entradas tenga, es 
    //decir es vector que se ajusta pues, 
    std::iota(indices.begin(), indices.end(), 0);
        // Llena un rango de elementos con valores secuenciales que comienzan desde un valor inicial.
        
         // llena ese sort orden a los ordenes que tenga
    std::sort(indices.begin(), indices.end(),
              [&pt](size_t a, size_t b) { return pt[a] > pt[b]; });



    // Aqui los ordena acorde al pT a los indices pues, ejemplo fueran 3 jets 30  10  50  los indices ordenados quedaria 3 1 2  
    
    // Aqui defino la funcion que ordena, toma los vectores originales con temp, siendo vec el vector a redordenar, sea pT, Eta, Phi, Mass
    
    auto reorder = [&indices](auto& vec) {
        auto temp = vec;
        for (size_t i = 0; i < indices.size(); ++i) {
            vec[i] = temp[indices[i]];       // aqui es acorde a los indices guardados en por el sort, es crea un nuevo vector que reorganiza los datos acorde como se organizaron
        }
    };

    // aqui es donde se le indica quienes son las que le va a hacer eso pues, que es como se definieron las entradas
    
    reorder(pt);
    reorder(eta);
    reorder(phi);
    reorder(mass);
}

// ===================================================
// Es aqui donde se tiene lo importante, que es bucle que llamara todas las funciones anteriores y las usa, aqui ya esta aplicado para que todos se reduzcan sin importar su tamaño de jets
// ===================================================

// Aqui se definen como const porque los va a tomar directamente del RDataframe por ende no lo va a modificar

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> processJets(const ROOT::VecOps::RVec<float>& jet_pt,
                                                          const ROOT::VecOps::RVec<float>& jet_eta,
                                                          const ROOT::VecOps::RVec<float>& jet_phi,
                                                          const ROOT::VecOps::RVec<float>& jet_mass) {
                                                          
    ROOT::VecOps::RVec<float> pt = jet_pt;
    ROOT::VecOps::RVec<float> eta = jet_eta;
    ROOT::VecOps::RVec<float> phi = jet_phi;
    ROOT::VecOps::RVec<float> mass = jet_mass;

    

    while (pt.size() > 3) {

        auto [idx0, idx1] = findClosestJets(phi, eta, pt, mass);
        
        TLorentzVector jet0, jet1, closest0, closest1;
        jet0.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
        jet1.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
        closest0.SetPtEtaPhiM(pt[idx0], eta[idx0], phi[idx0], mass[idx0]);
        closest1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], mass[idx1]);

        TLorentzVector sumajet0 = jet0 + closest0;
        TLorentzVector sumajet1 = jet1 + closest1;

        pt[0] = sumajet0.Pt();
        eta[0] = sumajet0.Eta();
        phi[0] = sumajet0.Phi();
        mass[0] = sumajet0.M();

        pt[1] = sumajet1.Pt();
        eta[1] = sumajet1.Eta();
        phi[1] = sumajet1.Phi();
        mass[1] = sumajet1.M();

        int idx_max = std::max(idx0, idx1);
        int idx_min = std::min(idx0, idx1);

        pt.erase(pt.begin() + idx_max);
        eta.erase(eta.begin() + idx_max);
        phi.erase(phi.begin() + idx_max);
        mass.erase(mass.begin() + idx_max);

        pt.erase(pt.begin() + idx_min);
        eta.erase(eta.begin() + idx_min);
        phi.erase(phi.begin() + idx_min);
        mass.erase(mass.begin() + idx_min);

        sortJetsByPt(pt, eta, phi, mass);


  
    }

    // ================= Cuando quedan 3 jets ====================
    if (pt.size() == 3) {
        
        TLorentzVector jet0, jet1, jet2;
        jet0.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
        jet1.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
        jet2.SetPtEtaPhiM(pt[2], eta[2], phi[2], mass[2]);

        double dr0 = jet2.DeltaR(jet0); 
        double dr1 = jet2.DeltaR(jet1);
        
        if (dr0 < dr1) {
            jet0 += jet2;
            pt[0] = jet0.Pt();
            eta[0] = jet0.Eta();
            phi[0] = jet0.Phi();
            mass[0] = jet0.M();

            pt.erase(pt.begin() + 2);
            eta.erase(eta.begin() + 2);
            phi.erase(phi.begin() + 2);
            mass.erase(mass.begin() + 2);
        } else {
            jet1 += jet2;
            pt[1] = jet1.Pt();
            eta[1] = jet1.Eta();
            phi[1] = jet1.Phi();
            mass[1] = jet1.M();

            pt.erase(pt.begin() + 2);
            eta.erase(eta.begin() + 2);
            phi.erase(phi.begin() + 2);
            mass.erase(mass.begin() + 2);
        }

        sortJetsByPt(pt, eta, phi, mass);


    }

    sortJetsByPt(pt, eta, phi, mass);

    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (size_t i = 0; i < pt.size(); ++i) {
        result.push_back({pt[i], eta[i], phi[i], mass[i]});
    }

    return result;
}


// ===================================================
// Todos los eventos
// ===================================================
void codigofull() {
    std::string input_dir = "data/";
    std::string output_dir = "plots/";
    std::string root_file = input_dir + "BL_M100_Pt3.root";

    std::string base_name = std::filesystem::path(root_file).stem().string(); // "M30_PT1"  

    ROOT::RDataFrame df("Delphes", root_file);

    auto df_processed = df
        .Define("processed_jets", processJets, {"Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
        .Define("NjetsFinales", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) {
            return static_cast<UInt_t>(jets.size());
        }, {"processed_jets"});

    auto df_twojets = df_processed.Filter("NjetsFinales >= 2")
        .Define("Masa_invariante", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) -> float {
            TLorentzVector jet1, jet2;
            jet1.SetPtEtaPhiM(jets[0][0], jets[0][1], jets[0][2], jets[0][3]);
            jet2.SetPtEtaPhiM(jets[1][0], jets[1][1], jets[1][2], jets[1][3]);
            return static_cast<float>((jet1 + jet2).M());
        }, {"processed_jets"})
        .Define("pt_lider",  [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[0][0]; }, {"processed_jets"})
        .Define("pt_sublider",  [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[1][0]; }, {"processed_jets"})
        .Define("eta_lider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[0][1]; }, {"processed_jets"})
        .Define("eta_sublider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[1][1]; }, {"processed_jets"})
        .Define("phi_lider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[0][2]; }, {"processed_jets"})
        .Define("phi_sublider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) { return jets[1][2]; }, {"processed_jets"});

    auto h_masa = df_twojets.Histo1D({"h_masa", "Masa invariante de 2 jets;Masa invariante [GeV];Eventos", 600, 0.0, 200}, "Masa_invariante");

    // Histogramas diferenciando líder y sublíder
    auto h_pt_lider     = df_twojets.Histo1D({"h_pt_lider", "pT del jet lider; pT [GeV];Eventos", 60, 0, 100}, "pt_lider");
    auto h_pt_sublider  = df_twojets.Histo1D({"h_pt_sublider", "pT del jet sublider; pT [GeV];Eventos", 60, 0, 100}, "pt_sublider");

    auto h_eta_lider    = df_twojets.Histo1D({"h_eta_lider", "eta del jet lider; eta;Eventos", 60, -6, 6}, "eta_lider");
    auto h_eta_sublider = df_twojets.Histo1D({"h_eta_sublider", "eta del jet sublider; eta;Eventos", 60, -6, 6}, "eta_sublider");

    auto h_phi_lider    = df_twojets.Histo1D({"h_phi_lider", "phi del jet lider; phi [rad];Eventos", 60, -4, 4}, "phi_lider");
    auto h_phi_sublider = df_twojets.Histo1D({"h_phi_sublider", "phi del jet sublider; phi [rad];Eventos", 60, -4, 4}, "phi_sublider");

    // Dibujar todos los histogramas
    TCanvas c1("c1", "Masa invariante", 800, 600);
    h_masa->Draw();
    c1.SaveAs((output_dir + base_name + "_h_masa_invariante.png").c_str());

    TCanvas c2("c2", "pT jets", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.12);
    h_pt_lider->SetLineColor(kRed);
    h_pt_sublider->SetLineColor(kBlue);
    double maxY = std::max(h_pt_lider->GetMaximum(), h_pt_sublider->GetMaximum());
    h_pt_lider->SetMaximum(1.2 * maxY);
    h_pt_lider->Draw();
    h_pt_sublider->Draw("SAME");
    auto legend = new TLegend(0.40, 0.65, 0.60, 0.88);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_pt_lider.GetPtr(), "pT del jet lider", "l");
    legend->AddEntry(h_pt_sublider.GetPtr(), "pT del jet sublider", "l");
    legend->Draw();
    c2.SaveAs((output_dir + base_name + "_h_pt_jets.png").c_str());

    TCanvas c3("c3", "Eta jets", 800, 600);
    h_eta_lider->SetLineColor(kRed);
    h_eta_sublider->SetLineColor(kBlue);
    h_eta_lider->Draw();
    h_eta_sublider->Draw("SAME");
    c3.BuildLegend();
    c3.SaveAs((output_dir + base_name + "_h_eta_jets.png").c_str());

    TCanvas c4("c4", "Phi jets", 800, 600);
    h_phi_lider->SetLineColor(kRed);
    h_phi_sublider->SetLineColor(kBlue);
    h_phi_lider->Draw();
    h_phi_sublider->Draw("SAME");
    c4.BuildLegend();
    c4.SaveAs((output_dir + base_name + "_h_phi_jets.png").c_str());

    // Histograma de número de jets por evento
    auto h_num_jets = df_processed.Histo1D(
        {"h_num_jets", "Numero de jets por evento post-reduccion;# de jets;Eventos", 10, -0.5, 100},
        "NjetsFinales"
    );
    TCanvas c5("c5", "Jets por evento", 800, 600);
    h_num_jets->SetFillColor(kAzure + 2);
    h_num_jets->SetLineColor(kBlack);
    h_num_jets->Draw("HIST");
    c5.SaveAs((output_dir + base_name + "_h_num_jets_por_evento.png").c_str());

    // pT combinado
    auto df_all_pt = df_twojets.Define("pt_lider_sublider", [](float pt1, float pt2) {
        return ROOT::VecOps::RVec<float>{pt1, pt2};
    }, {"pt_lider", "pt_sublider"});
    auto h_pt_all = df_all_pt.Histo1D({"h_pt_all", "Distribucion combinada de pT; pT [GeV]; Frecuencia", 60, 0, 100}, "pt_lider_sublider");
    TCanvas c6("c6", "pT combinado", 800, 600);
    h_pt_all->SetFillColor(kGreen + 2);
    h_pt_all->SetLineColor(kBlack);
    h_pt_all->Draw("HIST");
    c6.SaveAs((output_dir + base_name + "_h_pt_combinado.png").c_str());

    // eta combinado
    auto df_all_eta = df_all_pt.Define("eta_lider_sublider", [](float eta1, float eta2) {
        return ROOT::VecOps::RVec<float>{eta1, eta2};
    }, {"eta_lider", "eta_sublider"});
    auto h_eta_all = df_all_eta.Histo1D({"h_eta_all", "Distribucion combinada de #eta; #eta; Frecuencia", 60, -6, 6}, "eta_lider_sublider");
    TCanvas c7("c7", "#eta combinado", 800, 600);
    h_eta_all->SetFillColor(kOrange + 1);
    h_eta_all->SetLineColor(kBlack);
    h_eta_all->Draw("HIST");
    c7.SaveAs((output_dir + base_name + "_h_eta_combinado.png").c_str());

    // phi combinado
    auto df_all_phi = df_all_eta.Define("phi_lider_sublider", [](float phi1, float phi2) {
        return ROOT::VecOps::RVec<float>{phi1, phi2};
    }, {"phi_lider", "phi_sublider"});
    auto h_phi_all = df_all_phi.Histo1D({"h_phi_all", "Distribucion combinada de #phi; #phi [rad]; Frecuencia", 60, -4, 4}, "phi_lider_sublider");
    TCanvas c8("c8", "#phi combinado", 800, 600);
    h_phi_all->SetFillColor(kMagenta + 2);
    h_phi_all->SetLineColor(kBlack);
    h_phi_all->Draw("HIST");
    c8.SaveAs((output_dir + base_name + "_h_phi_combinado.png").c_str());

    // masa líder
    auto df_mass_lider = df_twojets.Define("mass_lider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) {
        return jets[0][3];
    }, {"processed_jets"});
    auto h_mass_lider = df_mass_lider.Histo1D({"h_mass_lider", "Masa del jet lider;Masa [GeV];Eventos", 60, 0.0, 100}, "mass_lider");
    TCanvas c9("c9", "Masa del jet lider", 800, 600);
    h_mass_lider->SetFillColor(kViolet + 1);
    h_mass_lider->SetLineColor(kBlack);
    h_mass_lider->Draw("HIST");
    c9.SaveAs((output_dir + base_name + "_h_masa_jet_lider.png").c_str());

    // masa sublíder
    auto df_mass_sublider = df_twojets.Define("mass_sublider", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) {
        return jets[1][3];
    }, {"processed_jets"});
    auto h_mass_sublider = df_mass_sublider.Histo1D({"h_mass_sublider", "Masa del jet sublider;Masa [GeV];Eventos", 60, 0.0, 50}, "mass_sublider");
    TCanvas c10("c10", "Masa del jet sublider", 800, 600);
    h_mass_sublider->SetFillColor(kOrange + 7);
    h_mass_sublider->SetLineColor(kBlack);
    h_mass_sublider->Draw("HIST");
    c10.SaveAs((output_dir + base_name + "_h_masa_jet_sublider.png").c_str());

    auto df_sistema = df_processed
        .Define("Masa_sistema", [](const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>& jets) -> float {
            TLorentzVector sistema;
            for (auto &j : jets) {
                TLorentzVector temp;
                temp.SetPtEtaPhiM(j[0], j[1], j[2], j[3]);
                sistema += temp;
            }
            return static_cast<float>(sistema.M());
        }, {"processed_jets"});

    // Histograma de masa del sistema completo
    auto h_masa_sistema = df_sistema.Histo1D(
        {"h_masa_sistema", "Masa invariante del sistema; Masa [GeV]; Eventos", 50, 0, 300},
        "Masa_sistema"
    );

    // Dibujar
    TCanvas c11("c11", "Masa del sistema", 800, 600);
    h_masa_sistema->SetLineColor(kBlack);
    h_masa_sistema->SetFillColor(kTeal + 2);
    h_masa_sistema->Draw("HIST");
    c11.SaveAs((output_dir + base_name + "_h_masa_sistema.png").c_str());


    
}
