#Function convert.clone.to.branch
convert.clone.to.branch <- function(t, branch.lens = NULL,
                                    merged.tree.node.annotation='sample.with.nonzero.cell.frac.ci'){
  # symbols used for each branch
  syms = c(seq(1,9), unlist(strsplit('abcdefghijklmnopqrstuvwxyz', '')))
  t = t[!is.na(t$parent) & !is.na(t$excluded) & !t$excluded,]
  t$branches = NA
  t$blengths = NA
  rownames(t) = t$lab
  grow.tree <- function(t, lab, parent.symbol=''){
    #print(paste0('---', lab))
    if (t[lab,'parent'] == '-1'){
      t[lab,'branches'] = 'Y'
    }
    children = !is.na(t$parent) & t$parent == lab
    if (any(children)){
      children.labs = t$lab[children]
      num.children = length(children.labs)
      if (num.children == 1){# && parent.symbol != ''){
        # one child, grow straight branch (no left right)
        children.syms = '0'
      }else{
        children.syms = syms[1:num.children]
      }
      children.syms = paste0(parent.symbol, children.syms)
      t$branches[children] = children.syms
      #print(children.labs)
      #print(children.syms)
      for (i in 1:length(children.labs)){
        t = grow.tree(t, children.labs[i], children.syms[i])
      }
    }
    return(t)
  }
  tg = grow.tree(t, t$lab[!is.na(t$parent) & t$parent == '-1'])
  if (merged.tree.node.annotation=='sample.with.nonzero.cell.frac.ci'){
    # remove ci info from sample annotation
    tg$samples.with.nonzero.cell.frac = gsub(',+$', '',
                                             gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.nonzero.cell.frac.ci))
    #tg$samples.with.nonzero.cell.frac = gsub(',+$', '', gsub('\u00B0[^,]+(,|$)', '',
    #     gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.cell.frac.ci)))
  }else{
    cat(paste0('WARN: merged.tree.node.annotation = ',
               merged.tree.node.annotation, ' not supported! No node annotation made.\n'))
  }
  if (is.null(branch.lens)){
    tg$blengths = 5
  }else{
    tg$blengths = branch.lens[tg$lab]
  }
  # color founding clone of met with diff. border
  #tg$node.border.color = ifelse(
  #    grepl('*', gsub('*P', '', tg[[merged.tree.node.annotation]], fixed=TRUE), fixed=T),
  #    'red', 'black')
  tg$node.border.color = 'black'
  tg$node.border.width = 1
  tg$branch.border.color = 'white'
  tg$branch.border.linetype = 'solid'
  tg$branch.border.width = 0.5
  
  return(tg)
  
}