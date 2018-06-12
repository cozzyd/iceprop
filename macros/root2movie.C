
struct movie_opts
{
  bool log = true; 
  double min = 1e-8; 
  double max = 1; 
  const char * tmpdir = ".movie"; 
  int width = 1920; 
  int height = 1080; 
  int max_entry = -1; 
  bool line_at_zero = true; 


}; 

TH2 * toFeet(TH2* h)
{

  const double meters_to_ft = 3.28084; 
  TH2 * hfeet = new TH2D(TString::Format("%s_ft",h->GetName()), h->GetTitle(), 
                         h->GetNbinsX(), h->GetXaxis()->GetXmin() * meters_to_ft, h->GetXaxis()->GetXmax()*meters_to_ft, 
                         h->GetNbinsY(), h->GetYaxis()->GetXmin() * meters_to_ft, h->GetYaxis()->GetXmax()*meters_to_ft); 

  TString xaxis(h->GetXaxis()->GetTitle()); 
  xaxis.ReplaceAll("r (m)","d (ft)"); 
  TString yaxis(h->GetYaxis()->GetTitle()); 
  yaxis.ReplaceAll("(m)","(ft)"); 
  hfeet->GetXaxis()->SetTitle(xaxis); 
  hfeet->GetYaxis()->SetTitle(yaxis); 


  for (int i = 1; i < h->GetNbinsX(); i++) 
  {
    for (int j = 1; j < h->GetNbinsY(); j++) 
    {
      hfeet->SetBinContent(i,j, h->GetBinContent(i,j)); 
    }
  }

  return hfeet; 
}

void root2movie(const char* fname, const char * movie_name = "movie.mp4",
                const char * tree_name = "Ez_mag", 
                movie_opts opts = movie_opts() )
{

  bool batch = gROOT->IsBatch(); 
  gROOT->SetBatch(true); 

  TFile f(fname); 
  TTree * tree = (TTree*) f.Get(tree_name); 


  TH2 * hist = 0; 
  tree->SetBranchAddress("hist",&hist); 


  TCanvas * c = new TCanvas("c1","canvas", opts.width, opts.height); 
  c->SetWindowSize(opts.width,opts.height); 

  if (opts.log) c->SetLogz(); 

  TString str; 
  str.Form("mkdir -p %s", opts.tmpdir); 
  system(str.Data()); 

  int maxi = opts.max_entry;
  if (maxi < 0) maxi = tree->GetEntries();

  TGraph line(2); 

  for (int i = 0; i < maxi; i++)
  {
    tree->GetEntry(i); 
    TH2 * ft = toFeet(hist); 

    if (opts.min)
    {
      ft->SetMinimum(opts.min); 
    }
    if (opts.max) 
    {
      ft->SetMaximum(opts.max); 
    }

    ft->Draw("colz"); 
    ft->SetStats(false); 

    if (opts.line_at_zero) 
    {
      line.SetPoint(0,ft->GetXaxis()->GetXmin(),0); 
      line.SetPoint(1,ft->GetXaxis()->GetXmax(),0); 
      line.Draw("lsame"); 
    }


    str.Form("%s/frame_%06d.png", opts.tmpdir, i); 
    c->SaveAs(str.Data()); 
    delete ft; 
  }

  str.Form("ffmpeg -y -i %s/frame_%%06d.png -c:v libx264 -pix_fmt yuv420p %s", opts.tmpdir, movie_name); 
  system(str.Data()); 
  str.Form("rm -rf %s", opts.tmpdir); 
  system(str.Data()); 

  delete hist; 
  gROOT->SetBatch(batch); 
}
