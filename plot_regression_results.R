# Reads gauss.RData, illustrations.RData and sensitivity.Rdata.
# Produces five PNG figures and four LaTeX tables printed to the console.
library(ggplot2)
library(xtable)
# Run from the directory containing gauss.RData, illustrations.RData and sensitivity.Rdata.
# Required object names and summary columns are checked explicitly.
folder <- getwd()
read_saved_file <- function(filename,object_name) {
  path <- file.path(folder,filename)
  if(!file.exists(path))stop('Missing file: ',path)
  saved <- new.env(parent=emptyenv())
  object_names <- load(path,envir=saved)
  if(!identical(object_names,object_name))
  stop(filename,' must contain exactly the object ',object_name,'. Rerun the current simulation script.')
  object <- saved[[object_name]]
  if(!is.list(object))stop('Expected a list named ',object_name,' in ',filename)
  object
}

# Missing components mean that their study was deliberately disabled in the run.
get_summary <- function(name) {
  sensitivity_components <- c('bandwidth','blocks','constants','coordinate')
  if(name %in% sensitivity_components) {
    component <- sensitivity[[name]]
  }else if(name=='comparison') {
    component <- list(summary=gauss$summary)
  }else if(name=='illustrations') {
    component <- illustrations$distributions
  }else {
    component <- illustrations[[name]]
  }
  if(is.null(component)) {
    message('Study not run: ',name,'; skipping its output.')
    return(NULL)
  }
  summary <- component$summary
  required <- c('n','d','epsilon','null_reps','alt_reps','type_I','type_II',
    'power','total_error','mcse')
  if(!is.data.frame(summary) || !all(required %in% names(summary)))
  stop('Malformed summary for ',name)
  if(any(!is.finite(as.matrix(summary[required]))))
  stop('Non-finite values in summary for ',name)
  if(any(abs(summary$total_error-summary$type_I-summary$type_II)>1e-10))
  stop('Inconsistent error probabilities for ',name)
  summary
}
check_unique <- function(s,keys)stopifnot(!anyDuplicated(s[keys]))
check_fixed <- function(s,keys)stopifnot(all(vapply(s[keys],function(x)length(unique(x))==1L,logical(1))))

# Widths are in cm; change illustration_width to adjust both illustration plots.
comparison_width <- 18
illustration_width <- 10
plot_error <- function(s,group,colours,labels,filename,width=illustration_width,height=7.2,
                 legend_order=names(colours)) {
  if(is.null(s)||!nrow(s))return(invisible(NULL))
  check_unique(s,c('n',group))
  stopifnot(all(as.character(s[[group]]) %in% names(colours)))
  s$curve <- factor(s[[group]],levels=names(colours))
  # The experiment script stores the Monte Carlo SE for total error.
  s$se <- s$mcse
  s$lower <- pmax(0,s$total_error-s$se)
  s$upper <- pmin(2,s$total_error+s$se)
  shapes <- setNames(rep(c(16,15,17,18,8),length.out=length(colours)),names(colours))
  p <- ggplot(s,aes(n,total_error,colour=curve,shape=curve,group=curve))+
  geom_line(linewidth=.65)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=40,linewidth=.35,show.legend=FALSE)+
  geom_point(size=2)+
  scale_colour_manual(values=colours,breaks=legend_order,
    labels=labels[match(legend_order,names(colours))])+
  scale_shape_manual(values=shapes,breaks=legend_order,
    labels=labels[match(legend_order,names(colours))])+
  scale_x_continuous(breaks=pretty(s$n,n=if(width>10)5 else 3))+
  scale_y_continuous(limits=c(0,max(1,1.05*max(s$upper))),expand=expansion(mult=c(.015,.03)))+
  labs(x=expression('Sample size '*n),y='Sum of type-I and -II error',colour=NULL,shape=NULL)+
  theme_classic(base_size=if(width>10)11 else 9.5)+
  theme(legend.position='right',legend.text=element_text(size=rel(1)),
    axis.text=element_text(colour='#333333'),
    panel.grid.major.y=element_line(colour='#E9E9E9'),plot.margin=margin(5,5,5,5))
  ggsave(file.path(folder,filename),plot=p,width=width,height=height,units='cm',
    dpi=320,device=grDevices::png,bg='white')
}
# All sensitivity tables report summed errors to three decimal places.
# xtable prints captions before labels; the top rule precedes grouped headings.
print_table <- function(s,rows,column,row_labels,heading,title,label) {
  if(is.null(s)||!nrow(s))return(invisible(NULL))
  check_unique(s,c(rows,column))
  check_fixed(s,c('n','d','epsilon','null_reps','alt_reps'))
  columns <- sort(unique(s[[column]]))
  if(column=='method')columns <- c('ToT-F','PB-LRT')[c('ToT-F','PB-LRT') %in% columns]
  s$total_error <- s$type_I+s$type_II
  wide <- reshape(s[c(rows,column,'total_error')],idvar=rows,timevar=column,direction='wide')
  wide <- wide[do.call(order,wide[rows]),c(rows,paste0('total_error.',columns)),drop=FALSE]
  nr <- length(rows)
  nc <- length(columns)
  names(wide) <- c(row_labels,if(is.numeric(s[[column]]))
    paste0('$',format(columns,trim=TRUE,scientific=FALSE),'$')else as.character(columns))
  header <- paste0('\\toprule\n',strrep(' & ',nr),
    '\\multicolumn{',nc,'}{c}{',heading,'} \\\\\n',
    '\\cmidrule(lr){',nr+1,'-',nr+nc,'}\n')
  caption <- sprintf(paste0('%s $n=%g$, $d=%g$, $\\epsilon=%g$. ',
    'Entries give estimated type-I + type-II error, using %d null and %d alternative repetitions per cell.'),
    title,s$n[1],s$d[1],s$epsilon[1],s$null_reps[1],s$alt_reps[1])
  print(xtable(wide,caption=caption,label=label,
    digits=c(0,ifelse(rows %in% c('C_sep','blocks'),0,2),rep(3,nc))),
    type='latex',file='',print.results=TRUE,include.rownames=FALSE,
    sanitize.colnames.function=identity,booktabs=TRUE,caption.placement='top',
    hline.after=c(0,nrow(wide)),add.to.row=list(pos=list(-1L),command=header),
    NA.string='--',comment=FALSE)
  cat('\n')
}

# 1. run_gaussian.R -> gauss.RData: Gaussian competitor comparison.
gauss <- read_saved_file('gauss.RData','gauss')
method_colours <- c('Transformed mean'='#2166AC','dpHSIC'='#B2182B',
  'ToT-F'='#238B45','PB-LRT'='#984EA3','Coordinate test'='#E08214')
plot_error(get_summary('comparison'),'method',method_colours,names(method_colours),
  'comparison.png',width=comparison_width,height=10.5,
  legend_order=c('Transformed mean','Coordinate test','dpHSIC','ToT-F','PB-LRT'))

# 2. run_illustrations.R -> illustrations.RData: all four illustration studies.
illustrations <- read_saved_file('illustrations.RData','illustrations')
s <- get_summary('illustrations')
if(!is.null(s))for(model in intersect(c('Cauchy','Exponential'),unique(s$model))) {
  g <- s[s$model==model,,drop=FALSE]
  budgets <- sort(unique(g$epsilon))
  colours <- setNames(rep(unname(method_colours),length.out=length(budgets)),as.character(budgets))
  labels <- as.expression(lapply(budgets,function(ep)bquote(epsilon == .(ep))))
  plot_error(g,'epsilon',colours,labels,paste0(tolower(model),'.png'))
}

# Four sparsity curves, with exactly the same tuning at every sparsity and n.
s <- get_summary('coordinate_sparsity')
if(!is.null(s)) {
  check_fixed(s,c('d','epsilon','eta'))
  sparsities <- sort(unique(s$sparsity))
  colours <- setNames(c('#2166AC','#B2182B','#238B45','#984EA3')[seq_along(sparsities)],as.character(sparsities))
  labels <- as.expression(lapply(sparsities,function(k)bquote(paste("\u2016", beta, "\u2016")[0] == .(k))))
  plot_error(s,'sparsity',colours,labels,'coordinate_sparsity.png',
    width=comparison_width,height=10.5)
}

# Response-only test: dimension on a logarithmic x axis.
s <- get_summary('response_dimension')
if(!is.null(s)) {
  check_unique(s,c('d','epsilon'))
  check_fixed(s,c('n'))
  s$curve <- factor(s$epsilon,levels=sort(unique(s$epsilon)))
  s$lower <- pmax(0,s$total_error-s$mcse)
  s$upper <- pmin(2,s$total_error+s$mcse)
  budgets <- sort(unique(s$epsilon))
  labels <- as.expression(lapply(budgets,function(ep)bquote(epsilon == .(ep))))
  p <- ggplot(s,aes(d,total_error,colour=curve,shape=curve,group=curve))+
    geom_line(linewidth=.65)+geom_point(size=2)+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.025,linewidth=.35)+
    scale_colour_manual(values=c('#2166AC','#B2182B','#238B45','#984EA3'),labels=labels)+
    scale_shape_manual(values=c(16,15,17,18),labels=labels)+
    scale_x_log10(breaks=sort(unique(s$d)))+
    labs(x='Dimension',y='Sum of type-I and -II error',colour=NULL,shape=NULL)+
    theme_classic(base_size=11)+theme(legend.position='right',
      legend.text=element_text(size=rel(1)))
  ggsave(file.path(folder,'dimension_stability.png'),plot=p,width=18,height=10.5,
         units='cm',dpi=320,bg='white')
}

# 3. run_sensitivity.R -> sensitivity.Rdata: all sensitivity tables.
sensitivity <- read_saved_file('sensitivity.Rdata','sensitivity')
print_table(get_summary('bandwidth'),'x_multiplier','sigma_Y',
  '$\\sigma_X/\\sqrt{d}$','$\\sigma_Y$',
  'Sensitivity of dpHSIC to kernel bandwidths.','tab:dpHSIC_sensitivity')
print_table(get_summary('blocks'),'blocks','method','Blocks','Method',
  'Sensitivity to the number of aggregation blocks.','tab:aggregate_blocks_sensitivity')

# One constants table: (C_sep,C_clip) rows; C_tol columns under each model heading.
s <- get_summary('constants')
if(!is.null(s)) {
  C_tol_values <- sort(unique(s$C_tol))
  k <- length(C_tol_values)
  models <- c('Gaussian','Exponential')
  if(!all(models %in% s$model)) {
    message('Constants table requires both completed models; skipping until both are available.')
  }else {
    for(model in models)check_fixed(s[s$model==model,],c('n','d','epsilon','eta','null_reps','alt_reps'))
    make_wide <- function(model,suffix) {
      g <- s[s$model==model,c('C_sep','C_clip','C_tol','type_I','type_II')]
      g$total_error <- g$type_I+g$type_II
      g <- g[c('C_sep','C_clip','C_tol','total_error')]
      check_unique(g,c('C_sep','C_clip','C_tol'))
      w <- reshape(g,idvar=c('C_sep','C_clip'),timevar='C_tol',direction='wide')
      value_columns <- paste0('total_error.',C_tol_values)
      for(name in setdiff(value_columns,names(w)))w[[name]] <- NA_real_
      w <- w[c('C_sep','C_clip',value_columns)]
      names(w)[-(1:2)] <- paste0('error',seq_len(k),suffix)
      w
    }
    tab <- merge(make_wide('Gaussian','_G'),make_wide('Exponential','_E'),by=c('C_sep','C_clip'),all=TRUE)
    tab <- tab[order(tab$C_sep,tab$C_clip),c('C_sep','C_clip',paste0('error',seq_len(k),'_G'),paste0('error',seq_len(k),'_E'))]
    names(tab) <- c('$C_{\\mathrm{sep}}$','$C_{\\mathrm{clip}}$',rep(paste0('$C_{\\mathrm{tol}}=',C_tol_values,'$'),2))
    header <- paste0('\\toprule\n & & \\multicolumn{',k,'}{c}{Gaussian linear regression}',
      ' & \\multicolumn{',k,'}{c}{Exponential regression} \\\\\n',
      '\\cmidrule(lr){3-',2+k,'}\\cmidrule(lr){',3+k,'-',2+2*k,'}\n')
    settings_text <- vapply(models,function(model) {
        g <- s[s$model==model,]
        sprintf('%s: $n=%g$, $d=%g$, $\\epsilon=%g$, $\\eta=%g$; %d null and %d alternative repetitions.',
        model,g$n[1],g$d[1],g$epsilon[1],g$eta[1],g$null_reps[1],g$alt_reps[1])
    },character(1))
    caption <- paste('Sensitivity to $C_{\\mathrm{sep}}$, $C_{\\mathrm{clip}}$ and $C_{\\mathrm{tol}}$. Entries give estimated type-I + type-II error.',
      paste(settings_text,collapse=' '),'Here $\\alpha=\\eta/(3J)$ in $l_s$.')
    print(xtable(tab,caption=caption,label='tab:sensitivity-constants',digits=c(0,0,2,rep(3,2*k))),
      type='latex',file='',print.results=TRUE,include.rownames=FALSE,
      sanitize.colnames.function=identity,caption.placement='top',comment=FALSE,
      NA.string='--',hline.after=NULL,add.to.row=list(pos=list(-1,0,nrow(tab)),
      command=c(header,'\\midrule\n','\\bottomrule\n')))
  }
}

# Coordinate sensitivity: sum of type-I and type-II errors, no bold entries.
# Both axes come from the saved grid (default: 0.1, 0.25, 1, 5, 10).
print_coordinate_table <- function(s) {
  if(is.null(s)||!nrow(s))return(invisible(NULL))
  check_fixed(s,c('n','d','epsilon','sparsity','null_reps','alt_reps'))
  check_unique(s,c('c_lambda','c_B'))
  cl <- sort(unique(s$c_lambda))
  cb <- sort(unique(s$c_B))
  caption <- sprintf(paste0(
    'Sensitivity of the coordinate test: $n=%g$, $d=%g$, ',
    '$\\epsilon=%g$, $\\|\\beta\\|_0=%g$. The working bounds are ',
    '$C_X=c_B\\sqrt{8/3}$ and $\\overline\\lambda=c_\\lambda$. ',
    'Entries give the sum of type-I and type-II errors from %d null and ',
    '%d alternative repetitions per cell.'),
    s$n[1],s$d[1],s$epsilon[1],s$sparsity[1],s$null_reps[1],s$alt_reps[1])
  row_end <- paste0("\\\\","\n")
  cat('\\begin{table}[ht]\n\\centering\n\\caption{',caption,
      '}\n\\label{tab:coordinate-sensitivity}\n',sep='')
  cat('\\begin{tabular}{',strrep('r',length(cb)+1),'}\n\\toprule\n',sep='')
  cat(' & \\multicolumn{',length(cb),'}{c}{$c_B$} ',row_end,
      '\\cmidrule(lr){2-',length(cb)+1,'}\n',sep='')
  cat('$c_\\lambda$ & ',
      paste0('$',format(cb,trim=TRUE,scientific=FALSE),'$',collapse=' & '),
      ' ',row_end,'\\midrule\n',sep='')
  for(a in cl) {
    entries <- vapply(cb,function(b) {
      cell <- s[s$c_lambda==a & s$c_B==b,,drop=FALSE]
      if(!nrow(cell))return('--')
      sprintf('%.3f',cell$type_I+cell$type_II)
    },character(1))
    cat(format(a,trim=TRUE,scientific=FALSE),' & ',
        paste(entries,collapse=' & '),' ',row_end,sep='')
  }
  cat('\\bottomrule\n\\end{tabular}\n\\end{table}\n')
}
print_coordinate_table(get_summary('coordinate'))
