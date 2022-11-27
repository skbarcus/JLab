void Test_Multiline_If()
{
  Double_t x = 10.;
  if(x==10)
    {
      cout<<Form("Single line worked!")<<endl;
    }

  if(
     x==5
     ||
     x==10
    )
    {
      cout<<Form("Multi line worked!")<<endl;
    }
}
