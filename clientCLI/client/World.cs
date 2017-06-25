﻿using System;
using System.Runtime.InteropServices;
using System.IO;
using System.Diagnostics;
using WorldCLR;
using System.Windows.Forms;

namespace client
{
    class World
    {
        private bool isInit = false;
        private IntPtr _world = IntPtr.Zero;
        [DllImport("kernel32.dll")]
        public static extern IntPtr LoadLibrary(string dllToLoad);

        [DllImport("kernel32.dll")]
        public static extern bool FreeLibrary(IntPtr hModule);

        public void init()
        {
            isInit = true;
            if (!Directory.Exists("players"))
                Directory.CreateDirectory("players");

            foreach (var curFile in Directory.GetFiles("players"))
                File.Delete(curFile);
        }

        public void addPlayers(params int[] submits)
        {
        }

        public void addBots(params string[] names)
        {
        }

        private static RunWorld rw;

        public void start()
        {
            try
            {
                /*if (_world != IntPtr.Zero)
                    destroyWorld();
                _world = LoadLibrary("World.dll");
                */
                try
                {
                    rw = new RunWorld();
                }
                catch(Exception ex)
                {
                    MessageBox.Show(ex.Message);
                    MessageBox.Show(ex.StackTrace);
                    
                }
            }
            catch (Exception ex)
            {
                Logger.LogMessage(ex);
            }
        }

        private void destroyWorld()
        {
            freePlayers();
            var a = FreeLibrary(_world);
        }

        public void freePlayers()
        {
            foreach (ProcessModule mod in Process.GetCurrentProcess().Modules)
            {
                try
                {
                    if ((new DirectoryInfo(mod.FileName)).Parent.ToString() == "players")
                        FreeLibrary(mod.BaseAddress);
                }
                catch (Exception ex)
                { }
            }
        }
    }
}
