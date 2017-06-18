using System;
using System.Runtime.InteropServices;
using System.IO;
using System.Diagnostics;

namespace client
{
    class World
    {
        public Server server;
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
            if (!isInit) init();
            foreach (var submit in submits)
                server.downloadDll(submit);
        }

        public void addBots(params string[] names)
        {
            if (!isInit) init();
            foreach (var bot in names)
                server.downloadBotDll(bot);
        }

        public void start()
        {
            try
            {
                if (_world != IntPtr.Zero)
                    destroyWorld();
                _world = LoadLibrary("include\\world.dll");
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
