setlocal tabstop=4 softtabstop=4 shiftwidth=4 expandtab

"settings :: Nvim-R plugin
nmap <LocalLeader>\ :RSend 
vmap <Space> <Plug>RDSendSelection
let R_assign=2

"R output is highlighted with current colorscheme
let g:rout_follow_colorscheme = 1

"R commands in R output are highlighted
let g:Rout_more_colors = 1

"let r_syntax_folding = 1
setlocal foldmethod=marker
