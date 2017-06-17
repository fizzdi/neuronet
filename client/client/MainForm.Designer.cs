namespace client
{
    partial class MainForm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(MainForm));
            this.tabcontrol = new System.Windows.Forms.TabControl();
            this.tab_send = new System.Windows.Forms.TabPage();
            this.b_clear_send = new System.Windows.Forms.Button();
            this.b_send = new System.Windows.Forms.Button();
            this.l_sourceFile = new System.Windows.Forms.LinkLabel();
            this.b_downloadFile = new System.Windows.Forms.Button();
            this.rtf_sourse = new System.Windows.Forms.RichTextBox();
            this.tab_history = new System.Windows.Forms.TabPage();
            this.b_reload_history = new System.Windows.Forms.Button();
            this.dgv_history = new System.Windows.Forms.DataGridView();
            this.col_timeOfSubmit = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.col_numberOfSubmit = new System.Windows.Forms.DataGridViewLinkColumn();
            this.col_submitState = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.col_showButton = new System.Windows.Forms.DataGridViewButtonColumn();
            this.tab_local_run = new System.Windows.Forms.TabPage();
            this.b_localRun = new System.Windows.Forms.Button();
            this.gb_local_selectedStrategy = new System.Windows.Forms.GroupBox();
            this.lb_local_selectedStrategy = new System.Windows.Forms.ListBox();
            this.b_local_left = new System.Windows.Forms.Button();
            this.b_local_right = new System.Windows.Forms.Button();
            this.gb_local_bots = new System.Windows.Forms.GroupBox();
            this.lb_local_bots = new System.Windows.Forms.ListBox();
            this.gb_local_myStrategy = new System.Windows.Forms.GroupBox();
            this.lb_local_myStrategy = new System.Windows.Forms.ListBox();
            this.tab_global_run = new System.Windows.Forms.TabPage();
            this.b_globalRun = new System.Windows.Forms.Button();
            this.gb_global_selectedStrategy = new System.Windows.Forms.GroupBox();
            this.lb_global_selectedStrategy = new System.Windows.Forms.ListBox();
            this.b_global_left = new System.Windows.Forms.Button();
            this.b_global_right = new System.Windows.Forms.Button();
            this.gb_global_otherStrategy = new System.Windows.Forms.GroupBox();
            this.lb_global_otherStrategy = new System.Windows.Forms.ListBox();
            this.gb_global_myStrategy = new System.Windows.Forms.GroupBox();
            this.lb_global_myStrategy = new System.Windows.Forms.ListBox();
            this.tab_settings = new System.Windows.Forms.TabPage();
            this.b_save_settings = new System.Windows.Forms.Button();
            this.tb_GroupName = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.ofd_sourseFile = new System.Windows.Forms.OpenFileDialog();
            this.button1 = new System.Windows.Forms.Button();
            this.tabcontrol.SuspendLayout();
            this.tab_send.SuspendLayout();
            this.tab_history.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.dgv_history)).BeginInit();
            this.tab_local_run.SuspendLayout();
            this.gb_local_selectedStrategy.SuspendLayout();
            this.gb_local_bots.SuspendLayout();
            this.gb_local_myStrategy.SuspendLayout();
            this.tab_global_run.SuspendLayout();
            this.gb_global_selectedStrategy.SuspendLayout();
            this.gb_global_otherStrategy.SuspendLayout();
            this.gb_global_myStrategy.SuspendLayout();
            this.tab_settings.SuspendLayout();
            this.SuspendLayout();
            // 
            // tabcontrol
            // 
            this.tabcontrol.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.tabcontrol.Controls.Add(this.tab_send);
            this.tabcontrol.Controls.Add(this.tab_history);
            this.tabcontrol.Controls.Add(this.tab_local_run);
            this.tabcontrol.Controls.Add(this.tab_global_run);
            this.tabcontrol.Controls.Add(this.tab_settings);
            this.tabcontrol.Location = new System.Drawing.Point(12, 12);
            this.tabcontrol.Name = "tabcontrol";
            this.tabcontrol.SelectedIndex = 0;
            this.tabcontrol.Size = new System.Drawing.Size(485, 337);
            this.tabcontrol.TabIndex = 0;
            this.tabcontrol.SelectedIndexChanged += new System.EventHandler(this.tabcontrol_SelectedIndexChanged);
            // 
            // tab_send
            // 
            this.tab_send.Controls.Add(this.button1);
            this.tab_send.Controls.Add(this.b_clear_send);
            this.tab_send.Controls.Add(this.b_send);
            this.tab_send.Controls.Add(this.l_sourceFile);
            this.tab_send.Controls.Add(this.b_downloadFile);
            this.tab_send.Controls.Add(this.rtf_sourse);
            this.tab_send.Location = new System.Drawing.Point(4, 22);
            this.tab_send.Name = "tab_send";
            this.tab_send.Padding = new System.Windows.Forms.Padding(3);
            this.tab_send.Size = new System.Drawing.Size(477, 311);
            this.tab_send.TabIndex = 0;
            this.tab_send.Text = "Отправка";
            this.tab_send.UseVisualStyleBackColor = true;
            // 
            // b_clear_send
            // 
            this.b_clear_send.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_clear_send.Location = new System.Drawing.Point(377, 36);
            this.b_clear_send.Name = "b_clear_send";
            this.b_clear_send.Size = new System.Drawing.Size(94, 23);
            this.b_clear_send.TabIndex = 4;
            this.b_clear_send.Text = "Очистить";
            this.b_clear_send.UseVisualStyleBackColor = true;
            this.b_clear_send.Click += new System.EventHandler(this.b_clear_send_Click);
            // 
            // b_send
            // 
            this.b_send.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_send.Location = new System.Drawing.Point(377, 7);
            this.b_send.Name = "b_send";
            this.b_send.Size = new System.Drawing.Size(94, 23);
            this.b_send.TabIndex = 3;
            this.b_send.Text = "Отправить";
            this.b_send.UseVisualStyleBackColor = true;
            this.b_send.Click += new System.EventHandler(this.b_send_Click);
            // 
            // l_sourceFile
            // 
            this.l_sourceFile.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.l_sourceFile.AutoSize = true;
            this.l_sourceFile.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.l_sourceFile.LinkColor = System.Drawing.SystemColors.ControlText;
            this.l_sourceFile.Location = new System.Drawing.Point(122, 282);
            this.l_sourceFile.Name = "l_sourceFile";
            this.l_sourceFile.Size = new System.Drawing.Size(92, 13);
            this.l_sourceFile.TabIndex = 2;
            this.l_sourceFile.TabStop = true;
            this.l_sourceFile.Text = "Файл не выбран";
            this.l_sourceFile.Click += new System.EventHandler(this.b_downloadFile_Click);
            // 
            // b_downloadFile
            // 
            this.b_downloadFile.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.b_downloadFile.Location = new System.Drawing.Point(7, 277);
            this.b_downloadFile.Name = "b_downloadFile";
            this.b_downloadFile.Size = new System.Drawing.Size(109, 23);
            this.b_downloadFile.TabIndex = 1;
            this.b_downloadFile.Text = "Выберите файл";
            this.b_downloadFile.UseVisualStyleBackColor = true;
            this.b_downloadFile.Click += new System.EventHandler(this.b_downloadFile_Click);
            // 
            // rtf_sourse
            // 
            this.rtf_sourse.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.rtf_sourse.Location = new System.Drawing.Point(7, 7);
            this.rtf_sourse.Name = "rtf_sourse";
            this.rtf_sourse.Size = new System.Drawing.Size(363, 264);
            this.rtf_sourse.TabIndex = 0;
            this.rtf_sourse.Text = "";
            // 
            // tab_history
            // 
            this.tab_history.Controls.Add(this.b_reload_history);
            this.tab_history.Controls.Add(this.dgv_history);
            this.tab_history.Location = new System.Drawing.Point(4, 22);
            this.tab_history.Name = "tab_history";
            this.tab_history.Padding = new System.Windows.Forms.Padding(3);
            this.tab_history.Size = new System.Drawing.Size(477, 311);
            this.tab_history.TabIndex = 1;
            this.tab_history.Text = "История";
            this.tab_history.UseVisualStyleBackColor = true;
            // 
            // b_reload_history
            // 
            this.b_reload_history.Location = new System.Drawing.Point(7, 7);
            this.b_reload_history.Name = "b_reload_history";
            this.b_reload_history.Size = new System.Drawing.Size(75, 23);
            this.b_reload_history.TabIndex = 1;
            this.b_reload_history.Text = "Обновить";
            this.b_reload_history.UseVisualStyleBackColor = true;
            this.b_reload_history.Click += new System.EventHandler(this.b_reload_history_Click);
            // 
            // dgv_history
            // 
            this.dgv_history.AllowUserToAddRows = false;
            this.dgv_history.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.dgv_history.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dgv_history.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.col_timeOfSubmit,
            this.col_numberOfSubmit,
            this.col_submitState,
            this.col_showButton});
            this.dgv_history.EditMode = System.Windows.Forms.DataGridViewEditMode.EditProgrammatically;
            this.dgv_history.Location = new System.Drawing.Point(7, 36);
            this.dgv_history.Name = "dgv_history";
            this.dgv_history.ReadOnly = true;
            this.dgv_history.RowHeadersVisible = false;
            this.dgv_history.RowHeadersWidthSizeMode = System.Windows.Forms.DataGridViewRowHeadersWidthSizeMode.DisableResizing;
            this.dgv_history.SelectionMode = System.Windows.Forms.DataGridViewSelectionMode.FullRowSelect;
            this.dgv_history.ShowCellErrors = false;
            this.dgv_history.ShowCellToolTips = false;
            this.dgv_history.ShowEditingIcon = false;
            this.dgv_history.ShowRowErrors = false;
            this.dgv_history.Size = new System.Drawing.Size(464, 269);
            this.dgv_history.TabIndex = 0;
            this.dgv_history.CellContentClick += new System.Windows.Forms.DataGridViewCellEventHandler(this.dgv_history_CellContentClick);
            this.dgv_history.Resize += new System.EventHandler(this.dgv_history_Resize);
            // 
            // col_timeOfSubmit
            // 
            this.col_timeOfSubmit.HeaderText = "Время отправки";
            this.col_timeOfSubmit.MinimumWidth = 60;
            this.col_timeOfSubmit.Name = "col_timeOfSubmit";
            this.col_timeOfSubmit.ReadOnly = true;
            this.col_timeOfSubmit.Resizable = System.Windows.Forms.DataGridViewTriState.True;
            this.col_timeOfSubmit.Width = 60;
            // 
            // col_numberOfSubmit
            // 
            this.col_numberOfSubmit.HeaderText = "Номер решения";
            this.col_numberOfSubmit.LinkColor = System.Drawing.SystemColors.ControlText;
            this.col_numberOfSubmit.MinimumWidth = 60;
            this.col_numberOfSubmit.Name = "col_numberOfSubmit";
            this.col_numberOfSubmit.ReadOnly = true;
            this.col_numberOfSubmit.Resizable = System.Windows.Forms.DataGridViewTriState.True;
            this.col_numberOfSubmit.SortMode = System.Windows.Forms.DataGridViewColumnSortMode.Automatic;
            this.col_numberOfSubmit.VisitedLinkColor = System.Drawing.SystemColors.ControlText;
            this.col_numberOfSubmit.Width = 60;
            // 
            // col_submitState
            // 
            this.col_submitState.HeaderText = "Статус";
            this.col_submitState.MinimumWidth = 100;
            this.col_submitState.Name = "col_submitState";
            this.col_submitState.ReadOnly = true;
            this.col_submitState.Width = 261;
            // 
            // col_showButton
            // 
            this.col_showButton.HeaderText = "Посмотреть";
            this.col_showButton.Name = "col_showButton";
            this.col_showButton.ReadOnly = true;
            this.col_showButton.Width = 80;
            // 
            // tab_local_run
            // 
            this.tab_local_run.Controls.Add(this.b_localRun);
            this.tab_local_run.Controls.Add(this.gb_local_selectedStrategy);
            this.tab_local_run.Controls.Add(this.b_local_left);
            this.tab_local_run.Controls.Add(this.b_local_right);
            this.tab_local_run.Controls.Add(this.gb_local_bots);
            this.tab_local_run.Controls.Add(this.gb_local_myStrategy);
            this.tab_local_run.Location = new System.Drawing.Point(4, 22);
            this.tab_local_run.Name = "tab_local_run";
            this.tab_local_run.Size = new System.Drawing.Size(477, 311);
            this.tab_local_run.TabIndex = 2;
            this.tab_local_run.Text = "Локальный запуск";
            this.tab_local_run.UseVisualStyleBackColor = true;
            this.tab_local_run.Resize += new System.EventHandler(this.tab_local_run_Resize);
            // 
            // b_localRun
            // 
            this.b_localRun.Anchor = System.Windows.Forms.AnchorStyles.Bottom;
            this.b_localRun.Enabled = true;
            this.b_localRun.Font = new System.Drawing.Font("Microsoft Sans Serif", 10F);
            this.b_localRun.Location = new System.Drawing.Point(207, 264);
            this.b_localRun.Name = "b_localRun";
            this.b_localRun.Size = new System.Drawing.Size(62, 44);
            this.b_localRun.TabIndex = 4;
            this.b_localRun.Text = "Запуск";
            this.b_localRun.UseVisualStyleBackColor = true;
            this.b_localRun.Click += new System.EventHandler(this.b_localRun_Click);
            // 
            // gb_local_selectedStrategy
            // 
            this.gb_local_selectedStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_local_selectedStrategy.Controls.Add(this.lb_local_selectedStrategy);
            this.gb_local_selectedStrategy.Location = new System.Drawing.Point(272, 3);
            this.gb_local_selectedStrategy.Name = "gb_local_selectedStrategy";
            this.gb_local_selectedStrategy.Size = new System.Drawing.Size(202, 305);
            this.gb_local_selectedStrategy.TabIndex = 4;
            this.gb_local_selectedStrategy.TabStop = false;
            this.gb_local_selectedStrategy.Text = "Выбранные стратегии:";
            // 
            // lb_local_selectedStrategy
            // 
            this.lb_local_selectedStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_local_selectedStrategy.FormattingEnabled = true;
            this.lb_local_selectedStrategy.Location = new System.Drawing.Point(6, 20);
            this.lb_local_selectedStrategy.Name = "lb_local_selectedStrategy";
            this.lb_local_selectedStrategy.Size = new System.Drawing.Size(190, 277);
            this.lb_local_selectedStrategy.TabIndex = 0;
            this.lb_local_selectedStrategy.SelectedIndexChanged += new System.EventHandler(this.lb_selectedStrategy_SelectedIndexChanged);
            // 
            // b_local_left
            // 
            this.b_local_left.Anchor = System.Windows.Forms.AnchorStyles.Bottom;
            this.b_local_left.Enabled = false;
            this.b_local_left.Font = new System.Drawing.Font("Microsoft Sans Serif", 14.25F);
            this.b_local_left.Location = new System.Drawing.Point(207, 179);
            this.b_local_left.Name = "b_local_left";
            this.b_local_left.Size = new System.Drawing.Size(62, 44);
            this.b_local_left.TabIndex = 3;
            this.b_local_left.Text = "<<";
            this.b_local_left.UseVisualStyleBackColor = true;
            this.b_local_left.Click += new System.EventHandler(this.b_localLeft_Click);
            // 
            // b_local_right
            // 
            this.b_local_right.Anchor = System.Windows.Forms.AnchorStyles.Top;
            this.b_local_right.Enabled = false;
            this.b_local_right.Font = new System.Drawing.Font("Microsoft Sans Serif", 14.25F);
            this.b_local_right.Location = new System.Drawing.Point(207, 76);
            this.b_local_right.Name = "b_local_right";
            this.b_local_right.Size = new System.Drawing.Size(62, 44);
            this.b_local_right.TabIndex = 2;
            this.b_local_right.Text = ">>";
            this.b_local_right.UseVisualStyleBackColor = true;
            this.b_local_right.Click += new System.EventHandler(this.b_localRight_Click);
            // 
            // gb_local_bots
            // 
            this.gb_local_bots.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_local_bots.Controls.Add(this.lb_local_bots);
            this.gb_local_bots.Location = new System.Drawing.Point(3, 159);
            this.gb_local_bots.Name = "gb_local_bots";
            this.gb_local_bots.Size = new System.Drawing.Size(201, 149);
            this.gb_local_bots.TabIndex = 1;
            this.gb_local_bots.TabStop = false;
            this.gb_local_bots.Text = "Боты:";
            // 
            // lb_local_bots
            // 
            this.lb_local_bots.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_local_bots.FormattingEnabled = true;
            this.lb_local_bots.Items.AddRange(new object[] {
            "I",
            "II",
            "III",
            "IV",
            "V"});
            this.lb_local_bots.Location = new System.Drawing.Point(7, 20);
            this.lb_local_bots.Name = "lb_local_bots";
            this.lb_local_bots.Size = new System.Drawing.Size(188, 121);
            this.lb_local_bots.Sorted = true;
            this.lb_local_bots.TabIndex = 0;
            this.lb_local_bots.Click += new System.EventHandler(this.lb_bots_Click);
            this.lb_local_bots.SelectedIndexChanged += new System.EventHandler(this.lb_bots_SelectedIndexChanged);
            // 
            // gb_local_myStrategy
            // 
            this.gb_local_myStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_local_myStrategy.Controls.Add(this.lb_local_myStrategy);
            this.gb_local_myStrategy.Location = new System.Drawing.Point(3, 3);
            this.gb_local_myStrategy.Name = "gb_local_myStrategy";
            this.gb_local_myStrategy.Size = new System.Drawing.Size(201, 150);
            this.gb_local_myStrategy.TabIndex = 0;
            this.gb_local_myStrategy.TabStop = false;
            this.gb_local_myStrategy.Text = "Ваши стратегии:";
            // 
            // lb_local_myStrategy
            // 
            this.lb_local_myStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_local_myStrategy.FormattingEnabled = true;
            this.lb_local_myStrategy.Items.AddRange(new object[] {
            "1",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9"});
            this.lb_local_myStrategy.Location = new System.Drawing.Point(7, 20);
            this.lb_local_myStrategy.Name = "lb_local_myStrategy";
            this.lb_local_myStrategy.Size = new System.Drawing.Size(188, 121);
            this.lb_local_myStrategy.Sorted = true;
            this.lb_local_myStrategy.TabIndex = 0;
            this.lb_local_myStrategy.Click += new System.EventHandler(this.lb_myStrategy_Click);
            this.lb_local_myStrategy.SelectedIndexChanged += new System.EventHandler(this.lb_myStrategy_SelectedIndexChanged);
            // 
            // tab_global_run
            // 
            this.tab_global_run.Controls.Add(this.b_globalRun);
            this.tab_global_run.Controls.Add(this.gb_global_selectedStrategy);
            this.tab_global_run.Controls.Add(this.b_global_left);
            this.tab_global_run.Controls.Add(this.b_global_right);
            this.tab_global_run.Controls.Add(this.gb_global_otherStrategy);
            this.tab_global_run.Controls.Add(this.gb_global_myStrategy);
            this.tab_global_run.Location = new System.Drawing.Point(4, 22);
            this.tab_global_run.Name = "tab_global_run";
            this.tab_global_run.Size = new System.Drawing.Size(477, 311);
            this.tab_global_run.TabIndex = 3;
            this.tab_global_run.Text = "Глобальный запуск";
            this.tab_global_run.UseVisualStyleBackColor = true;
            this.tab_global_run.Resize += new System.EventHandler(this.tab_global_run_Resize);
            // 
            // b_globalRun
            // 
            this.b_globalRun.Anchor = System.Windows.Forms.AnchorStyles.Bottom;
            this.b_globalRun.Enabled = false;
            this.b_globalRun.Font = new System.Drawing.Font("Microsoft Sans Serif", 10F);
            this.b_globalRun.Location = new System.Drawing.Point(207, 264);
            this.b_globalRun.Name = "b_globalRun";
            this.b_globalRun.Size = new System.Drawing.Size(62, 44);
            this.b_globalRun.TabIndex = 9;
            this.b_globalRun.Text = "Запуск";
            this.b_globalRun.UseVisualStyleBackColor = true;
            this.b_globalRun.Click += new System.EventHandler(this.b_globalRun_Click);
            // 
            // gb_global_selectedStrategy
            // 
            this.gb_global_selectedStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_global_selectedStrategy.Controls.Add(this.lb_global_selectedStrategy);
            this.gb_global_selectedStrategy.Location = new System.Drawing.Point(272, 3);
            this.gb_global_selectedStrategy.Name = "gb_global_selectedStrategy";
            this.gb_global_selectedStrategy.Size = new System.Drawing.Size(202, 305);
            this.gb_global_selectedStrategy.TabIndex = 10;
            this.gb_global_selectedStrategy.TabStop = false;
            this.gb_global_selectedStrategy.Text = "Выбранные стратегии:";
            // 
            // lb_global_selectedStrategy
            // 
            this.lb_global_selectedStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_global_selectedStrategy.FormattingEnabled = true;
            this.lb_global_selectedStrategy.Location = new System.Drawing.Point(6, 20);
            this.lb_global_selectedStrategy.Name = "lb_global_selectedStrategy";
            this.lb_global_selectedStrategy.Size = new System.Drawing.Size(190, 277);
            this.lb_global_selectedStrategy.TabIndex = 0;
            this.lb_global_selectedStrategy.SelectedIndexChanged += new System.EventHandler(this.lb_global_selectedStrategy_SelectedIndexChanged);
            // 
            // b_global_left
            // 
            this.b_global_left.Anchor = System.Windows.Forms.AnchorStyles.Bottom;
            this.b_global_left.Enabled = false;
            this.b_global_left.Font = new System.Drawing.Font("Microsoft Sans Serif", 14.25F);
            this.b_global_left.Location = new System.Drawing.Point(207, 179);
            this.b_global_left.Name = "b_global_left";
            this.b_global_left.Size = new System.Drawing.Size(62, 44);
            this.b_global_left.TabIndex = 8;
            this.b_global_left.Text = "<<";
            this.b_global_left.UseVisualStyleBackColor = true;
            this.b_global_left.Click += new System.EventHandler(this.b_global_left_Click);
            // 
            // b_global_right
            // 
            this.b_global_right.Anchor = System.Windows.Forms.AnchorStyles.Top;
            this.b_global_right.Enabled = false;
            this.b_global_right.Font = new System.Drawing.Font("Microsoft Sans Serif", 14.25F);
            this.b_global_right.Location = new System.Drawing.Point(207, 76);
            this.b_global_right.Name = "b_global_right";
            this.b_global_right.Size = new System.Drawing.Size(62, 44);
            this.b_global_right.TabIndex = 7;
            this.b_global_right.Text = ">>";
            this.b_global_right.UseVisualStyleBackColor = true;
            this.b_global_right.Click += new System.EventHandler(this.b_global_right_Click);
            // 
            // gb_global_otherStrategy
            // 
            this.gb_global_otherStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_global_otherStrategy.Controls.Add(this.lb_global_otherStrategy);
            this.gb_global_otherStrategy.Location = new System.Drawing.Point(3, 159);
            this.gb_global_otherStrategy.Name = "gb_global_otherStrategy";
            this.gb_global_otherStrategy.Size = new System.Drawing.Size(201, 149);
            this.gb_global_otherStrategy.TabIndex = 6;
            this.gb_global_otherStrategy.TabStop = false;
            this.gb_global_otherStrategy.Text = "Стратегии других участников:";
            // 
            // lb_global_otherStrategy
            // 
            this.lb_global_otherStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_global_otherStrategy.FormattingEnabled = true;
            this.lb_global_otherStrategy.Items.AddRange(new object[] {
            "I",
            "II",
            "III",
            "IV",
            "V"});
            this.lb_global_otherStrategy.Location = new System.Drawing.Point(7, 20);
            this.lb_global_otherStrategy.Name = "lb_global_otherStrategy";
            this.lb_global_otherStrategy.Size = new System.Drawing.Size(188, 121);
            this.lb_global_otherStrategy.Sorted = true;
            this.lb_global_otherStrategy.TabIndex = 0;
            this.lb_global_otherStrategy.Click += new System.EventHandler(this.lb_global_otherStrategy_Click);
            this.lb_global_otherStrategy.SelectedIndexChanged += new System.EventHandler(this.lb_global_otherStrategy_SelectedIndexChanged);
            // 
            // gb_global_myStrategy
            // 
            this.gb_global_myStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.gb_global_myStrategy.Controls.Add(this.lb_global_myStrategy);
            this.gb_global_myStrategy.Location = new System.Drawing.Point(3, 3);
            this.gb_global_myStrategy.Name = "gb_global_myStrategy";
            this.gb_global_myStrategy.Size = new System.Drawing.Size(201, 150);
            this.gb_global_myStrategy.TabIndex = 5;
            this.gb_global_myStrategy.TabStop = false;
            this.gb_global_myStrategy.Text = "Ваши стратегии:";
            // 
            // lb_global_myStrategy
            // 
            this.lb_global_myStrategy.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_global_myStrategy.FormattingEnabled = true;
            this.lb_global_myStrategy.Items.AddRange(new object[] {
            "1",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9"});
            this.lb_global_myStrategy.Location = new System.Drawing.Point(7, 20);
            this.lb_global_myStrategy.Name = "lb_global_myStrategy";
            this.lb_global_myStrategy.Size = new System.Drawing.Size(188, 121);
            this.lb_global_myStrategy.Sorted = true;
            this.lb_global_myStrategy.TabIndex = 0;
            this.lb_global_myStrategy.Click += new System.EventHandler(this.lb_global_myStrategy_Click);
            this.lb_global_myStrategy.SelectedIndexChanged += new System.EventHandler(this.lb_global_myStrategy_SelectedIndexChanged);
            // 
            // tab_settings
            // 
            this.tab_settings.Controls.Add(this.b_save_settings);
            this.tab_settings.Controls.Add(this.tb_GroupName);
            this.tab_settings.Controls.Add(this.label1);
            this.tab_settings.Location = new System.Drawing.Point(4, 22);
            this.tab_settings.Name = "tab_settings";
            this.tab_settings.Padding = new System.Windows.Forms.Padding(3);
            this.tab_settings.Size = new System.Drawing.Size(477, 311);
            this.tab_settings.TabIndex = 4;
            this.tab_settings.Text = "Настройки";
            this.tab_settings.UseVisualStyleBackColor = true;
            // 
            // b_save_settings
            // 
            this.b_save_settings.Location = new System.Drawing.Point(292, 59);
            this.b_save_settings.Name = "b_save_settings";
            this.b_save_settings.Size = new System.Drawing.Size(75, 23);
            this.b_save_settings.TabIndex = 2;
            this.b_save_settings.Text = "Применить";
            this.b_save_settings.UseVisualStyleBackColor = true;
            this.b_save_settings.Click += new System.EventHandler(this.b_save_settings_Click);
            // 
            // tb_GroupName
            // 
            this.tb_GroupName.Location = new System.Drawing.Point(136, 22);
            this.tb_GroupName.Name = "tb_GroupName";
            this.tb_GroupName.Size = new System.Drawing.Size(231, 20);
            this.tb_GroupName.TabIndex = 1;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(21, 25);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(109, 13);
            this.label1.TabIndex = 0;
            this.label1.Text = "Название команды:";
            // 
            // ofd_sourseFile
            // 
            this.ofd_sourseFile.Filter = "Файл кода C++|*.cpp";
            this.ofd_sourseFile.Title = "Выберите файл решения";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(377, 161);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 5;
            this.button1.Text = "button1";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(509, 361);
            this.Controls.Add(this.tabcontrol);
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.MinimumSize = new System.Drawing.Size(525, 399);
            this.Name = "MainForm";
            this.Text = "Client for Code Hungry Game Challenge ";
            this.Load += new System.EventHandler(this.MainForm_Load);
            this.tabcontrol.ResumeLayout(false);
            this.tab_send.ResumeLayout(false);
            this.tab_send.PerformLayout();
            this.tab_history.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.dgv_history)).EndInit();
            this.tab_local_run.ResumeLayout(false);
            this.gb_local_selectedStrategy.ResumeLayout(false);
            this.gb_local_bots.ResumeLayout(false);
            this.gb_local_myStrategy.ResumeLayout(false);
            this.tab_global_run.ResumeLayout(false);
            this.gb_global_selectedStrategy.ResumeLayout(false);
            this.gb_global_otherStrategy.ResumeLayout(false);
            this.gb_global_myStrategy.ResumeLayout(false);
            this.tab_settings.ResumeLayout(false);
            this.tab_settings.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.TabControl tabcontrol;
        private System.Windows.Forms.TabPage tab_send;
        private System.Windows.Forms.TabPage tab_history;
        private System.Windows.Forms.TabPage tab_local_run;
        private System.Windows.Forms.TabPage tab_global_run;
        private System.Windows.Forms.Button b_downloadFile;
        private System.Windows.Forms.RichTextBox rtf_sourse;
        private System.Windows.Forms.LinkLabel l_sourceFile;
        private System.Windows.Forms.OpenFileDialog ofd_sourseFile;
        private System.Windows.Forms.Button b_send;
        private System.Windows.Forms.Button b_clear_send;
        private System.Windows.Forms.Button b_reload_history;
        private System.Windows.Forms.DataGridView dgv_history;
        private System.Windows.Forms.DataGridViewTextBoxColumn col_timeOfSubmit;
        private System.Windows.Forms.DataGridViewLinkColumn col_numberOfSubmit;
        private System.Windows.Forms.DataGridViewTextBoxColumn col_submitState;
        private System.Windows.Forms.DataGridViewButtonColumn col_showButton;
        private System.Windows.Forms.GroupBox gb_local_selectedStrategy;
        private System.Windows.Forms.ListBox lb_local_selectedStrategy;
        private System.Windows.Forms.Button b_local_left;
        private System.Windows.Forms.Button b_local_right;
        private System.Windows.Forms.GroupBox gb_local_bots;
        private System.Windows.Forms.ListBox lb_local_bots;
        private System.Windows.Forms.GroupBox gb_local_myStrategy;
        private System.Windows.Forms.ListBox lb_local_myStrategy;
        private System.Windows.Forms.Button b_localRun;
        private System.Windows.Forms.Button b_globalRun;
        private System.Windows.Forms.GroupBox gb_global_selectedStrategy;
        private System.Windows.Forms.ListBox lb_global_selectedStrategy;
        private System.Windows.Forms.Button b_global_left;
        private System.Windows.Forms.Button b_global_right;
        private System.Windows.Forms.GroupBox gb_global_otherStrategy;
        private System.Windows.Forms.ListBox lb_global_otherStrategy;
        private System.Windows.Forms.GroupBox gb_global_myStrategy;
        private System.Windows.Forms.ListBox lb_global_myStrategy;
        private System.Windows.Forms.TabPage tab_settings;
        private System.Windows.Forms.TextBox tb_GroupName;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Button b_save_settings;
        private System.Windows.Forms.Button button1;
    }
}

