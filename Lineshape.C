#include <tpc2019/TpcPrototypeUnpacker.h>
R__LOAD_LIBRARY(libtpc2019.so)


double
SignalShape_PowerLawExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  double signal = par[0] * pow((x[0] - par[1]), par[2]) * exp(-(x[0] - par[1]) * par[3]);
  return pedestal + signal;
}

double
SignalShape_PowerLawDoubleExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;

  double t = (x[0] - par[1]);
  double tau = par[3] / par[2];
  double n = par[2];

//  double signal =                                                                                         //
//      par[0]                                                                                              //
//      * pow((x[0] - par[1]), par[2])                                                                      //
//      * ((1 / pow(par[3], par[2]) * exp(par[2])) * exp(-(x[0] - par[1]) * (par[2] / par[3]))  //   //
//        );
  double signal_core = pow(t/tau, n) * exp( - t/tau);
  double signal =signal_core* par[0] / (exp(-n) * pow(n, n));

  return pedestal + signal;
}

void Lineshape()
{
  gSystem->Load("libtpc2019");

  TFile *f_TpcPrototypeUnpacker = TFile::Open("/sphenix/user/jinhuang/TPC/fnal_June2019/SimPadPlaneIter7/tpc_beam_00000413-0000.evt_TpcPrototypeUnpacker.root");

  f_TpcPrototypeUnpacker->ls();

  TTree * eventT = nullptr;
  f_TpcPrototypeUnpacker->GetObject("eventT",eventT);

  eventT->Show(10);

  int event = 10;
  int verbosity = 2;

  TBranch *branch = eventT->GetBranch("Clusters");

  TClonesArray *m_IOClusters = new TClonesArray("TpcPrototypeUnpacker::ClusterData", 1000);
  branch->SetAddress(&m_IOClusters);

  eventT->GetEvent(event);

  cout << "Fetched " << m_IOClusters->GetEntries() << " clusters" << endl;

  m_IOClusters->Print();

  for (int i = 0; i < m_IOClusters->GetEntries(); i++)
  {
    TpcPrototypeUnpacker::ClusterData *cluster = dynamic_cast<TpcPrototypeUnpacker::ClusterData *>(m_IOClusters->At(i));

    if (cluster == nullptr)
    {
      cout << "missing cluster " << i << endl;
    }
    else
    {
      double peak = NAN;
      double peak_sample = NAN;
      double pedestal = NAN;

      vector<double> samples = cluster->sum_samples;

      static const int n_parameter = 5;

      // inital guesses
      int peakPos = 0.;

      //  assert(samples.size() == n_samples);
      const int n_samples = samples.size();

      TGraph gpulse(n_samples);
      for (int i = 0; i < n_samples; i++)
      {
        (gpulse.GetX())[i] = i;

        (gpulse.GetY())[i] = samples[i];
      }

      pedestal = gpulse.GetY()[0];  //(double) PEDESTAL;
      double peakval = pedestal;
      const double risetime = 1.5;

      for (int iSample = 0; iSample < n_samples - risetime * 3; iSample++)
      {
        if (abs(gpulse.GetY()[iSample] - pedestal) > abs(peakval - pedestal))
        {
          peakval = gpulse.GetY()[iSample];
          peakPos = iSample;
        }
      }
      peakval -= pedestal;

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - "
             << "pedestal = " << pedestal << ", "
             << "peakval = " << peakval << ", "
             << "peakPos = " << peakPos << endl;
      }

      // build default value
      struct default_values_t
      {
        default_values_t(double default_value, double min_value, double max_value)
          : def(default_value)
          , min(min_value)
          , max(max_value)
        {
        }
        double def;
        double min;
        double max;
      };

      vector<default_values_t> default_values(n_parameter, default_values_t(numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN()));

      default_values[0] = default_values_t(peakval * .7, peakval * -1.5, peakval * 1.5);
      default_values[1] = default_values_t(peakPos + risetime, peakPos - 7 * risetime, peakPos - risetime);
      default_values[2] = default_values_t(2., 1, 10.);
      default_values[3] = default_values_t(risetime, risetime * .2, risetime * 10);
      default_values[4] = default_values_t(pedestal, pedestal - abs(peakval), pedestal + abs(peakval));

      // fit function
      TF1 fits("f_SignalShape_PowerLawDoubleExp", SignalShape_PowerLawDoubleExp, 0., n_samples, n_parameter);
      fits.SetParNames("Amplitude", "Sample Start", "Power", "Peak Time 1", "Pedestal");

      for (int i = 0; i < n_parameter; ++i)
      {
          fits.SetParameter(i, default_values[i].def);

          if (default_values[i].min < default_values[i].max)
          {
            fits.SetParLimits(i, default_values[i].min, default_values[i].max);
          }
          else
          {
            fits.FixParameter(i, default_values[i].def);
          }

          if (verbosity)
          {
            cout << "SampleFit_PowerLawDoubleExp - parameter [" << i << "]: "
                 << "default value = " << default_values[i].def
                 << ", min value = " << default_values[i].min
                 << ", max value = " << default_values[i].max << endl;
          }
      }

      // gpulse.Fit("f_SignalShape_PowerLawDoubleExp", "V0");

      if (verbosity <= 1)
        gpulse.Fit(&fits, "QRN0W", "goff", 0., (double) n_samples);
      else
        gpulse.Fit(&fits, "RN0W+", "goff", 0., (double) n_samples);

      // store results
      pedestal = fits.GetParameter(4);

      const double peakpos1 = fits.GetParameter(3);
      const double peakpos2 = fits.GetParameter(6);
      double max_peakpos = fits.GetParameter(1) + (peakpos1 > peakpos2 ? peakpos1 : peakpos2);
      if (max_peakpos > n_samples - 1) max_peakpos = n_samples - 1;

      if (fits.GetParameter(0) > 0)
        peak_sample = fits.GetMaximumX(fits.GetParameter(1), max_peakpos);
      else
        peak_sample = fits.GetMinimumX(fits.GetParameter(1), max_peakpos);

      peak = fits.Eval(peak_sample) - pedestal;

      if (verbosity)
      {
        static int id = 0;
        ++id;

        string c_name(string("SampleFit_PowerLawDoubleExp_") + to_string(id));

        TCanvas *canvas = new TCanvas(
            c_name.c_str(), c_name.c_str());
        canvas->Update();

        TGraph *g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
        g_plot->SetTitle((string("ADC data and fit #") + to_string(id) + string(";Sample number [50 ns];ADC value")).c_str());

        fits.SetLineColor(kMagenta);
        fits.DrawClone("same");
        fits.Print();


        TGraph g_max(1);

        g_max.GetX()[0] = peak_sample;
        g_max.GetY()[0] = peak + pedestal;

        g_max.SetMarkerStyle(kFullCircle);
        g_max.SetMarkerSize(2);
        g_max.SetMarkerColor(kRed);

        g_max.DrawClone("p");

        TLegend * leg = new TLegend(0.5,0.5,0.95, 0.95);
        for (int p = 0; p<n_parameter; ++p)
        {
          leg->AddEntry("",Form("%s = %.3f", fits.GetParName(p), fits.GetParameter(p)), "");
        }
        leg->Draw();

        canvas->Update();
        canvas->Draw();
        canvas->Print((c_name + ".png").c_str());
      }

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - "
             << "peak_sample = " << peak_sample << ", "
             << "max_peakpos = " << max_peakpos << ", "
             << "Sample Start = " << fits.GetParameter(1) << ", "
             << "peak = " << peak << ", "
             << "pedestal = " << pedestal << endl;
      }
    }
  }
}
