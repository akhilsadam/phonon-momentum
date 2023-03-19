#!/bin/bash
#------------------------------------
#
# Simple install script for vampire
# on ubuntu linux > 18.04
#
#------------------------------------
echo "Installing vampire at /opt/vampire"
echo "This script must be run with sudo/root privilages"
#create proram directories
echo "Creating installation directory"
mkdir /opt/vampire
mkdir /opt/vampire/bin
echo "Copying files to installation directory"
cp vampire-serial /opt/vampire/bin/
echo "Installation complete"
echo "Please add the following line to your shell configuration file (usually ~/.bashrc):"
echo "export PATH=/opt/vampire/bin/:\$PATH"

