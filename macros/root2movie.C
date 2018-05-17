
struct movie_opts
{
  bool log = true; 
  double min = 1e-10; 
  double max = 1; 
  const char * tmpdir = ".movie"; 
  int width = 1920; 
  int height = 1080; 
  int max_entry = -1; 


}; 

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

  for (int i = 0; i < maxi; i++)
  {
    tree->GetEntry(i); 

    if (opts.min)
    {
      hist->SetMinimum(opts.min); 
    }
    if (opts.max) 
    {
      hist->SetMaximum(opts.max); 
    }

    hist->Draw("colz"); 
    str.Form("%s/frame_%06d.png", opts.tmpdir, i); 
    c->SaveAs(str.Data()); 
  }

  str.Form("ffmpeg -y -i %s/frame_%%06d.png -c:v libx264 -pix_fmt yuv420p %s", opts.tmpdir, movie_name); 
  system(str.Data()); 
  str.Form("rm -rf %s", opts.tmpdir); 
  system(str.Data()); 

  delete hist; 
  gROOT->SetBatch(batch); 
}
