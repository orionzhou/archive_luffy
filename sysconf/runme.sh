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

[[ -d .config/pip ]] || mkdir -p .config/pip
ln -sf $CONFIGDIR/pip.conf .config/pip/pip.conf

[[ -h .vim ]] && rm .vim
[[ -d .vim ]] || mkdir .vim
cd $CONFIGDIR/vim
vimfolders=$(ls -d */)
cd ~/.vim
for i in $vimfolders
do
  j=${i%%/}
  [[ -h $j ]] && rm $j
  ln -sf $CONFIGDIR/vim/$j $j
done
cd ..

[[ -d Library/Application\ Support/BBEdit/Language\ Modules ]] && cp $CONFIGDIR/R.plist Library/Application\ Support/BBEdit/Language\ Modules

echo " all done -> cheers "

