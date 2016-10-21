function export_pgf(data,path,name)
   %  export_pgf(data,path,name) 
   tbl_name = [path '/' name '.txt'];
   pgfplot_tbl(data,[],tbl_name,[],[]);
end 
