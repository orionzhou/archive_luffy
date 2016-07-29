#!/bin/bash

CONFIGDIR=`pwd`

cd

ln -sf $CONFIGDIR/bash_profile .bash_profile
ln -sf $CONFIGDIR/bashrc .bashrc
ln -sf $CONFIGDIR/git.conf .gitconfig
ln -sf $CONFIGDIR/Rprofile .Rprofile
ln -sf $CONFIGDIR/vimrc .vimrc
ln -sf $CONFIGDIR/dircolors .dircolors
ln -sf $CONFIGDIR/hg.conf .hg.conf

[[ -h .ssh ]] && rm .ssh
[[ -d .ssh ]] || mkdir .ssh
ln -sf $CONFIGDIR/ssh.conf .ssh/config
chmod 600 $CONFIGDIR/ssh.conf
[[ -h .ssh/authorized_keys ]] && rm .ssh/authorized_keys
cp -f $CONFIGDIR/authorized_keys .ssh/authorized_keys

[[ -h .vim ]] && rm .vim
[[ -d .vim ]] || mkdir .vim
[[ -h .vim/colors ]] && rm .vim/colors
[[ -h .vim/ftdetect ]] && rm .vim/ftdetect
[[ -h .vim/ftplugin ]] && rm .vim/ftplugin
ln -sf $CONFIGDIR/vim/colors .vim/colors
ln -sf $CONFIGDIR/vim/ftdetect .vim/ftdetect
ln -sf $CONFIGDIR/vim/ftplugin .vim/ftplugin
[[ -d Library/Application\ Support/BBEdit/Language\ Modules ]] && cp $CONFIGDIR/R.plist Library/Application\ Support/BBEdit/Language\ Modules

echo " all done -> cheers "

