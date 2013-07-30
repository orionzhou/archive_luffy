#!/bin/bash

CONFIGDIR=`pwd`

cd

ln -sf $CONFIGDIR/bashrc .bashrc
ln -sf $CONFIGDIR/git.conf .gitconfig
ln -sf $CONFIGDIR/Rprofile .Rprofile
ln -sf $CONFIGDIR/vimrc .vimrc
ln -sf $CONFIGDIR/dircolors .dircolors

[[ -h .ssh ]] && rm .ssh
[[ -d .ssh ]] || mkdir .ssh
ln -sf $CONFIGDIR/ssh.conf .ssh/config
chmod 600 $CONFIGDIR/ssh.conf

[[ -h .vim ]] && rm .vim
[[ -d .vim ]] || mkdir .vim
[[ -h .vim/colors ]] && rm .vim/colors
[[ -h .vim/ftdetect ]] && rm .vim/ftdetect
[[ -h .vim/ftplugin ]] && rm .vim/ftplugin
ln -sf $CONFIGDIR/vim/colors .vim/colors
ln -sf $CONFIGDIR/vim/ftdetect .vim/ftdetect
ln -sf $CONFIGDIR/vim/ftplugin .vim/ftplugin

echo " all done -> cheers "

