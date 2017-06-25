namespace client
{
    partial class Edit_player
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
            this.p_net = new System.Windows.Forms.Panel();
            this.groupBox2 = new System.Windows.Forms.GroupBox();
            this.label8 = new System.Windows.Forms.Label();
            this.nud_q_learn = new System.Windows.Forms.NumericUpDown();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.label1 = new System.Windows.Forms.Label();
            this.nud_rms_eps = new System.Windows.Forms.NumericUpDown();
            this.label2 = new System.Windows.Forms.Label();
            this.nud_rms_learnrate = new System.Windows.Forms.NumericUpDown();
            this.label3 = new System.Windows.Forms.Label();
            this.cb_train_type = new System.Windows.Forms.ComboBox();
            this.label4 = new System.Windows.Forms.Label();
            this.nud_rms_gamma = new System.Windows.Forms.NumericUpDown();
            this.groupBox4 = new System.Windows.Forms.GroupBox();
            this.label5 = new System.Windows.Forms.Label();
            this.label16 = new System.Windows.Forms.Label();
            this.nud_train_eps = new System.Windows.Forms.NumericUpDown();
            this.nud_max_test_count = new System.Windows.Forms.NumericUpDown();
            this.nud_sensors = new System.Windows.Forms.NumericUpDown();
            this.label17 = new System.Windows.Forms.Label();
            this.nud_hid_neurons = new System.Windows.Forms.NumericUpDown();
            this.label18 = new System.Windows.Forms.Label();
            this.nud_hid_layers = new System.Windows.Forms.NumericUpDown();
            this.label19 = new System.Windows.Forms.Label();
            this.nud_test_count = new System.Windows.Forms.NumericUpDown();
            this.label20 = new System.Windows.Forms.Label();
            this.nud_train_epoch = new System.Windows.Forms.NumericUpDown();
            this.label21 = new System.Windows.Forms.Label();
            this.nud_train_period = new System.Windows.Forms.NumericUpDown();
            this.label22 = new System.Windows.Forms.Label();
            this.nud_end_train = new System.Windows.Forms.NumericUpDown();
            this.label23 = new System.Windows.Forms.Label();
            this.label24 = new System.Windows.Forms.Label();
            this.cb_fun_act = new System.Windows.Forms.ComboBox();
            this.b_save = new System.Windows.Forms.Button();
            this.groupBox5 = new System.Windows.Forms.GroupBox();
            this.rb_solve = new System.Windows.Forms.RadioButton();
            this.rb_net = new System.Windows.Forms.RadioButton();
            this.p_solve = new System.Windows.Forms.Panel();
            this.l_filename = new System.Windows.Forms.Label();
            this.rtf_solve = new System.Windows.Forms.RichTextBox();
            this.b_open = new System.Windows.Forms.Button();
            this.b_clear = new System.Windows.Forms.Button();
            this.linkLabel1 = new System.Windows.Forms.LinkLabel();
            this.ofd_sourseFile = new System.Windows.Forms.OpenFileDialog();
            this.label25 = new System.Windows.Forms.Label();
            this.tb_name = new System.Windows.Forms.TextBox();
            this.p_net.SuspendLayout();
            this.groupBox2.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_q_learn)).BeginInit();
            this.groupBox1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_eps)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_learnrate)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_gamma)).BeginInit();
            this.groupBox4.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_eps)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_max_test_count)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_sensors)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_hid_neurons)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_hid_layers)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_test_count)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_epoch)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_period)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_end_train)).BeginInit();
            this.groupBox5.SuspendLayout();
            this.p_solve.SuspendLayout();
            this.SuspendLayout();
            // 
            // p_net
            // 
            this.p_net.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.p_net.Controls.Add(this.label25);
            this.p_net.Controls.Add(this.tb_name);
            this.p_net.Controls.Add(this.groupBox2);
            this.p_net.Controls.Add(this.groupBox1);
            this.p_net.Controls.Add(this.groupBox4);
            this.p_net.Location = new System.Drawing.Point(12, 83);
            this.p_net.Name = "p_net";
            this.p_net.Size = new System.Drawing.Size(743, 355);
            this.p_net.TabIndex = 1;
            // 
            // groupBox2
            // 
            this.groupBox2.Controls.Add(this.label8);
            this.groupBox2.Controls.Add(this.nud_q_learn);
            this.groupBox2.Location = new System.Drawing.Point(387, 196);
            this.groupBox2.Name = "groupBox2";
            this.groupBox2.Size = new System.Drawing.Size(353, 52);
            this.groupBox2.TabIndex = 3;
            this.groupBox2.TabStop = false;
            this.groupBox2.Text = "Параметры Q-Learning:";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(7, 24);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(140, 17);
            this.label8.TabIndex = 29;
            this.label8.Text = "Скорость обучения:";
            // 
            // nud_q_learn
            // 
            this.nud_q_learn.DecimalPlaces = 2;
            this.nud_q_learn.Increment = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_q_learn.Location = new System.Drawing.Point(185, 22);
            this.nud_q_learn.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_q_learn.Minimum = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_q_learn.Name = "nud_q_learn";
            this.nud_q_learn.Size = new System.Drawing.Size(155, 22);
            this.nud_q_learn.TabIndex = 13;
            this.nud_q_learn.Value = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.label1);
            this.groupBox1.Controls.Add(this.nud_rms_eps);
            this.groupBox1.Controls.Add(this.label2);
            this.groupBox1.Controls.Add(this.nud_rms_learnrate);
            this.groupBox1.Controls.Add(this.label3);
            this.groupBox1.Controls.Add(this.cb_train_type);
            this.groupBox1.Controls.Add(this.label4);
            this.groupBox1.Controls.Add(this.nud_rms_gamma);
            this.groupBox1.Location = new System.Drawing.Point(387, 42);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(353, 148);
            this.groupBox1.TabIndex = 2;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "Настройки алгоритма обучения:";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(7, 105);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(106, 17);
            this.label1.TabIndex = 33;
            this.label1.Text = "RMS точность:";
            // 
            // nud_rms_eps
            // 
            this.nud_rms_eps.DecimalPlaces = 4;
            this.nud_rms_eps.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_rms_eps.Location = new System.Drawing.Point(185, 103);
            this.nud_rms_eps.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_rms_eps.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            this.nud_rms_eps.Name = "nud_rms_eps";
            this.nud_rms_eps.Size = new System.Drawing.Size(155, 22);
            this.nud_rms_eps.TabIndex = 3;
            this.nud_rms_eps.Value = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(7, 80);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(172, 17);
            this.label2.TabIndex = 31;
            this.label2.Text = "RMS скорость обучения:";
            // 
            // nud_rms_learnrate
            // 
            this.nud_rms_learnrate.DecimalPlaces = 4;
            this.nud_rms_learnrate.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_rms_learnrate.Location = new System.Drawing.Point(185, 78);
            this.nud_rms_learnrate.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_rms_learnrate.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            this.nud_rms_learnrate.Name = "nud_rms_learnrate";
            this.nud_rms_learnrate.Size = new System.Drawing.Size(155, 22);
            this.nud_rms_learnrate.TabIndex = 2;
            this.nud_rms_learnrate.Value = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(7, 53);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(88, 17);
            this.label3.TabIndex = 29;
            this.label3.Text = "RMS Гамма:";
            // 
            // cb_train_type
            // 
            this.cb_train_type.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.cb_train_type.FormattingEnabled = true;
            this.cb_train_type.Items.AddRange(new object[] {
            "RMS",
            "RPROP"});
            this.cb_train_type.Location = new System.Drawing.Point(185, 21);
            this.cb_train_type.Name = "cb_train_type";
            this.cb_train_type.Size = new System.Drawing.Size(155, 24);
            this.cb_train_type.TabIndex = 0;
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(7, 25);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(104, 17);
            this.label4.TabIndex = 27;
            this.label4.Text = "Тип обучения:";
            // 
            // nud_rms_gamma
            // 
            this.nud_rms_gamma.DecimalPlaces = 2;
            this.nud_rms_gamma.Increment = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_rms_gamma.Location = new System.Drawing.Point(185, 51);
            this.nud_rms_gamma.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_rms_gamma.Minimum = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_rms_gamma.Name = "nud_rms_gamma";
            this.nud_rms_gamma.Size = new System.Drawing.Size(155, 22);
            this.nud_rms_gamma.TabIndex = 1;
            this.nud_rms_gamma.Value = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            // 
            // groupBox4
            // 
            this.groupBox4.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left)));
            this.groupBox4.Controls.Add(this.label5);
            this.groupBox4.Controls.Add(this.label16);
            this.groupBox4.Controls.Add(this.nud_train_eps);
            this.groupBox4.Controls.Add(this.nud_max_test_count);
            this.groupBox4.Controls.Add(this.nud_sensors);
            this.groupBox4.Controls.Add(this.label17);
            this.groupBox4.Controls.Add(this.nud_hid_neurons);
            this.groupBox4.Controls.Add(this.label18);
            this.groupBox4.Controls.Add(this.nud_hid_layers);
            this.groupBox4.Controls.Add(this.label19);
            this.groupBox4.Controls.Add(this.nud_test_count);
            this.groupBox4.Controls.Add(this.label20);
            this.groupBox4.Controls.Add(this.nud_train_epoch);
            this.groupBox4.Controls.Add(this.label21);
            this.groupBox4.Controls.Add(this.nud_train_period);
            this.groupBox4.Controls.Add(this.label22);
            this.groupBox4.Controls.Add(this.nud_end_train);
            this.groupBox4.Controls.Add(this.label23);
            this.groupBox4.Controls.Add(this.label24);
            this.groupBox4.Controls.Add(this.cb_fun_act);
            this.groupBox4.Location = new System.Drawing.Point(9, 43);
            this.groupBox4.Name = "groupBox4";
            this.groupBox4.Size = new System.Drawing.Size(372, 301);
            this.groupBox4.TabIndex = 1;
            this.groupBox4.TabStop = false;
            this.groupBox4.Text = "Настройки нейронной сети:";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(6, 18);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(155, 17);
            this.label5.TabIndex = 10;
            this.label5.Text = "Количество сенсоров:";
            // 
            // label16
            // 
            this.label16.AutoSize = true;
            this.label16.Location = new System.Drawing.Point(6, 130);
            this.label16.Name = "label16";
            this.label16.Size = new System.Drawing.Size(185, 17);
            this.label16.TabIndex = 26;
            this.label16.Text = "Размер тестовой выборки:";
            // 
            // nud_train_eps
            // 
            this.nud_train_eps.DecimalPlaces = 4;
            this.nud_train_eps.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_train_eps.Location = new System.Drawing.Point(208, 212);
            this.nud_train_eps.Maximum = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.nud_train_eps.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            this.nud_train_eps.Name = "nud_train_eps";
            this.nud_train_eps.Size = new System.Drawing.Size(155, 22);
            this.nud_train_eps.TabIndex = 7;
            this.nud_train_eps.Value = new decimal(new int[] {
            1,
            0,
            0,
            262144});
            // 
            // nud_max_test_count
            // 
            this.nud_max_test_count.Location = new System.Drawing.Point(208, 128);
            this.nud_max_test_count.Maximum = new decimal(new int[] {
            25000,
            0,
            0,
            0});
            this.nud_max_test_count.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_max_test_count.Name = "nud_max_test_count";
            this.nud_max_test_count.Size = new System.Drawing.Size(155, 22);
            this.nud_max_test_count.TabIndex = 4;
            this.nud_max_test_count.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // nud_sensors
            // 
            this.nud_sensors.Location = new System.Drawing.Point(208, 16);
            this.nud_sensors.Maximum = new decimal(new int[] {
            64,
            0,
            0,
            0});
            this.nud_sensors.Minimum = new decimal(new int[] {
            4,
            0,
            0,
            0});
            this.nud_sensors.Name = "nud_sensors";
            this.nud_sensors.Size = new System.Drawing.Size(155, 22);
            this.nud_sensors.TabIndex = 0;
            this.nud_sensors.Value = new decimal(new int[] {
            4,
            0,
            0,
            0});
            // 
            // label17
            // 
            this.label17.AutoSize = true;
            this.label17.Location = new System.Drawing.Point(6, 271);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(171, 17);
            this.label17.TabIndex = 24;
            this.label17.Text = "Тип функции активации:";
            // 
            // nud_hid_neurons
            // 
            this.nud_hid_neurons.Location = new System.Drawing.Point(208, 44);
            this.nud_hid_neurons.Maximum = new decimal(new int[] {
            512,
            0,
            0,
            0});
            this.nud_hid_neurons.Minimum = new decimal(new int[] {
            8,
            0,
            0,
            0});
            this.nud_hid_neurons.Name = "nud_hid_neurons";
            this.nud_hid_neurons.Size = new System.Drawing.Size(155, 22);
            this.nud_hid_neurons.TabIndex = 1;
            this.nud_hid_neurons.Value = new decimal(new int[] {
            8,
            0,
            0,
            0});
            // 
            // label18
            // 
            this.label18.AutoSize = true;
            this.label18.Location = new System.Drawing.Point(6, 242);
            this.label18.Name = "label18";
            this.label18.Size = new System.Drawing.Size(179, 17);
            this.label18.TabIndex = 22;
            this.label18.Text = "Ход окончания обучения:";
            // 
            // nud_hid_layers
            // 
            this.nud_hid_layers.Location = new System.Drawing.Point(208, 72);
            this.nud_hid_layers.Maximum = new decimal(new int[] {
            5,
            0,
            0,
            0});
            this.nud_hid_layers.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_hid_layers.Name = "nud_hid_layers";
            this.nud_hid_layers.Size = new System.Drawing.Size(155, 22);
            this.nud_hid_layers.TabIndex = 2;
            this.nud_hid_layers.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(6, 214);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(141, 17);
            this.label19.TabIndex = 21;
            this.label19.Text = "Точность обучения:";
            // 
            // nud_test_count
            // 
            this.nud_test_count.Location = new System.Drawing.Point(208, 100);
            this.nud_test_count.Maximum = new decimal(new int[] {
            200,
            0,
            0,
            0});
            this.nud_test_count.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_test_count.Name = "nud_test_count";
            this.nud_test_count.Size = new System.Drawing.Size(155, 22);
            this.nud_test_count.TabIndex = 3;
            this.nud_test_count.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label20
            // 
            this.label20.AutoSize = true;
            this.label20.Location = new System.Drawing.Point(6, 186);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(129, 17);
            this.label20.TabIndex = 20;
            this.label20.Text = "Период обучения:";
            // 
            // nud_train_epoch
            // 
            this.nud_train_epoch.Location = new System.Drawing.Point(208, 156);
            this.nud_train_epoch.Maximum = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.nud_train_epoch.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_train_epoch.Name = "nud_train_epoch";
            this.nud_train_epoch.Size = new System.Drawing.Size(155, 22);
            this.nud_train_epoch.TabIndex = 5;
            this.nud_train_epoch.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label21
            // 
            this.label21.AutoSize = true;
            this.label21.Location = new System.Drawing.Point(6, 158);
            this.label21.Name = "label21";
            this.label21.Size = new System.Drawing.Size(110, 17);
            this.label21.TabIndex = 19;
            this.label21.Text = "Эпох обучения:";
            // 
            // nud_train_period
            // 
            this.nud_train_period.Location = new System.Drawing.Point(208, 184);
            this.nud_train_period.Maximum = new decimal(new int[] {
            5000,
            0,
            0,
            0});
            this.nud_train_period.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_train_period.Name = "nud_train_period";
            this.nud_train_period.Size = new System.Drawing.Size(155, 22);
            this.nud_train_period.TabIndex = 6;
            this.nud_train_period.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label22
            // 
            this.label22.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left)));
            this.label22.AutoSize = true;
            this.label22.Location = new System.Drawing.Point(6, 102);
            this.label22.Name = "label22";
            this.label22.Size = new System.Drawing.Size(176, 17);
            this.label22.TabIndex = 18;
            this.label22.Text = "Тестов в эпоху обучения:";
            // 
            // nud_end_train
            // 
            this.nud_end_train.Increment = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.nud_end_train.Location = new System.Drawing.Point(208, 240);
            this.nud_end_train.Maximum = new decimal(new int[] {
            500000,
            0,
            0,
            0});
            this.nud_end_train.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_end_train.Name = "nud_end_train";
            this.nud_end_train.Size = new System.Drawing.Size(155, 22);
            this.nud_end_train.TabIndex = 8;
            this.nud_end_train.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // label23
            // 
            this.label23.AutoSize = true;
            this.label23.Location = new System.Drawing.Point(6, 72);
            this.label23.Name = "label23";
            this.label23.Size = new System.Drawing.Size(111, 17);
            this.label23.TabIndex = 17;
            this.label23.Text = "Скрытых слоев:";
            // 
            // label24
            // 
            this.label24.AutoSize = true;
            this.label24.Location = new System.Drawing.Point(6, 46);
            this.label24.Name = "label24";
            this.label24.Size = new System.Drawing.Size(176, 17);
            this.label24.TabIndex = 11;
            this.label24.Text = "Нейронов скрытого слоя:";
            // 
            // cb_fun_act
            // 
            this.cb_fun_act.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.cb_fun_act.FormattingEnabled = true;
            this.cb_fun_act.Items.AddRange(new object[] {
            "Сигмойда",
            "Тангенс гиперболический"});
            this.cb_fun_act.Location = new System.Drawing.Point(208, 268);
            this.cb_fun_act.Name = "cb_fun_act";
            this.cb_fun_act.Size = new System.Drawing.Size(155, 24);
            this.cb_fun_act.TabIndex = 9;
            // 
            // b_save
            // 
            this.b_save.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.b_save.Location = new System.Drawing.Point(299, 444);
            this.b_save.Name = "b_save";
            this.b_save.Size = new System.Drawing.Size(130, 29);
            this.b_save.TabIndex = 2;
            this.b_save.Text = "Сохранить";
            this.b_save.UseVisualStyleBackColor = true;
            this.b_save.Click += new System.EventHandler(this.b_save_Click);
            // 
            // groupBox5
            // 
            this.groupBox5.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.groupBox5.Controls.Add(this.rb_solve);
            this.groupBox5.Controls.Add(this.rb_net);
            this.groupBox5.Location = new System.Drawing.Point(12, 24);
            this.groupBox5.Name = "groupBox5";
            this.groupBox5.Size = new System.Drawing.Size(756, 53);
            this.groupBox5.TabIndex = 0;
            this.groupBox5.TabStop = false;
            this.groupBox5.Text = "Вид стратегии игрока";
            // 
            // rb_solve
            // 
            this.rb_solve.AutoSize = true;
            this.rb_solve.Location = new System.Drawing.Point(157, 22);
            this.rb_solve.Name = "rb_solve";
            this.rb_solve.Size = new System.Drawing.Size(124, 21);
            this.rb_solve.TabIndex = 1;
            this.rb_solve.Text = "Свое решение";
            this.rb_solve.UseVisualStyleBackColor = true;
            this.rb_solve.CheckedChanged += new System.EventHandler(this.radioButton1_CheckedChanged);
            // 
            // rb_net
            // 
            this.rb_net.AutoSize = true;
            this.rb_net.Checked = true;
            this.rb_net.Location = new System.Drawing.Point(9, 22);
            this.rb_net.Name = "rb_net";
            this.rb_net.Size = new System.Drawing.Size(136, 21);
            this.rb_net.TabIndex = 0;
            this.rb_net.TabStop = true;
            this.rb_net.Text = "Нейронная сеть";
            this.rb_net.UseVisualStyleBackColor = true;
            // 
            // p_solve
            // 
            this.p_solve.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.p_solve.Controls.Add(this.l_filename);
            this.p_solve.Controls.Add(this.rtf_solve);
            this.p_solve.Controls.Add(this.b_open);
            this.p_solve.Controls.Add(this.b_clear);
            this.p_solve.Location = new System.Drawing.Point(12, 83);
            this.p_solve.Name = "p_solve";
            this.p_solve.Size = new System.Drawing.Size(743, 355);
            this.p_solve.TabIndex = 20;
            // 
            // l_filename
            // 
            this.l_filename.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.l_filename.AutoSize = true;
            this.l_filename.Location = new System.Drawing.Point(161, 329);
            this.l_filename.Name = "l_filename";
            this.l_filename.Size = new System.Drawing.Size(118, 17);
            this.l_filename.TabIndex = 5;
            this.l_filename.Text = "Файл не выбран";
            // 
            // rtf_solve
            // 
            this.rtf_solve.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.rtf_solve.Location = new System.Drawing.Point(4, 4);
            this.rtf_solve.Margin = new System.Windows.Forms.Padding(4);
            this.rtf_solve.Name = "rtf_solve";
            this.rtf_solve.Size = new System.Drawing.Size(593, 311);
            this.rtf_solve.TabIndex = 0;
            this.rtf_solve.Text = "";
            // 
            // b_open
            // 
            this.b_open.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.b_open.Location = new System.Drawing.Point(9, 323);
            this.b_open.Margin = new System.Windows.Forms.Padding(4);
            this.b_open.Name = "b_open";
            this.b_open.Size = new System.Drawing.Size(145, 28);
            this.b_open.TabIndex = 1;
            this.b_open.Text = "Выберите файл";
            this.b_open.UseVisualStyleBackColor = true;
            this.b_open.Click += new System.EventHandler(this.b_open_Click);
            // 
            // b_clear
            // 
            this.b_clear.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_clear.Location = new System.Drawing.Point(605, 13);
            this.b_clear.Margin = new System.Windows.Forms.Padding(4);
            this.b_clear.Name = "b_clear";
            this.b_clear.Size = new System.Drawing.Size(125, 28);
            this.b_clear.TabIndex = 4;
            this.b_clear.Text = "Очистить";
            this.b_clear.UseVisualStyleBackColor = true;
            this.b_clear.Click += new System.EventHandler(this.b_clear_Click);
            // 
            // linkLabel1
            // 
            this.linkLabel1.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.linkLabel1.AutoSize = true;
            this.linkLabel1.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.linkLabel1.LinkColor = System.Drawing.SystemColors.ControlText;
            this.linkLabel1.Location = new System.Drawing.Point(166, 340);
            this.linkLabel1.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.linkLabel1.Name = "linkLabel1";
            this.linkLabel1.Size = new System.Drawing.Size(0, 17);
            this.linkLabel1.TabIndex = 17;
            // 
            // ofd_sourseFile
            // 
            this.ofd_sourseFile.Filter = "Файл кода C++|*.cpp";
            this.ofd_sourseFile.Title = "Выберите файл решения";
            // 
            // label25
            // 
            this.label25.AutoSize = true;
            this.label25.Location = new System.Drawing.Point(18, 10);
            this.label25.Name = "label25";
            this.label25.Size = new System.Drawing.Size(152, 17);
            this.label25.TabIndex = 29;
            this.label25.Text = "Название персонажа:";
            // 
            // tb_name
            // 
            this.tb_name.Location = new System.Drawing.Point(182, 7);
            this.tb_name.MaxLength = 40;
            this.tb_name.Name = "tb_name";
            this.tb_name.Size = new System.Drawing.Size(272, 22);
            this.tb_name.TabIndex = 0;
            // 
            // Edit_player
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(780, 485);
            this.Controls.Add(this.b_save);
            this.Controls.Add(this.groupBox5);
            this.Controls.Add(this.linkLabel1);
            this.Controls.Add(this.p_net);
            this.Controls.Add(this.p_solve);
            this.Name = "Edit_player";
            this.Text = "Настройка персонажа";
            this.Load += new System.EventHandler(this.Edit_player_Load);
            this.p_net.ResumeLayout(false);
            this.p_net.PerformLayout();
            this.groupBox2.ResumeLayout(false);
            this.groupBox2.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_q_learn)).EndInit();
            this.groupBox1.ResumeLayout(false);
            this.groupBox1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_eps)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_learnrate)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_rms_gamma)).EndInit();
            this.groupBox4.ResumeLayout(false);
            this.groupBox4.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_eps)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_max_test_count)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_sensors)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_hid_neurons)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_hid_layers)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_test_count)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_epoch)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_train_period)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_end_train)).EndInit();
            this.groupBox5.ResumeLayout(false);
            this.groupBox5.PerformLayout();
            this.p_solve.ResumeLayout(false);
            this.p_solve.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Panel p_net;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.NumericUpDown nud_rms_eps;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.NumericUpDown nud_rms_learnrate;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.ComboBox cb_train_type;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.NumericUpDown nud_rms_gamma;
        private System.Windows.Forms.GroupBox groupBox4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label16;
        private System.Windows.Forms.NumericUpDown nud_train_eps;
        private System.Windows.Forms.NumericUpDown nud_max_test_count;
        private System.Windows.Forms.NumericUpDown nud_sensors;
        private System.Windows.Forms.Label label17;
        private System.Windows.Forms.NumericUpDown nud_hid_neurons;
        private System.Windows.Forms.Label label18;
        private System.Windows.Forms.NumericUpDown nud_hid_layers;
        private System.Windows.Forms.Label label19;
        private System.Windows.Forms.NumericUpDown nud_test_count;
        private System.Windows.Forms.Label label20;
        private System.Windows.Forms.NumericUpDown nud_train_epoch;
        private System.Windows.Forms.Label label21;
        private System.Windows.Forms.NumericUpDown nud_train_period;
        private System.Windows.Forms.Label label22;
        private System.Windows.Forms.NumericUpDown nud_end_train;
        private System.Windows.Forms.Label label23;
        private System.Windows.Forms.Label label24;
        private System.Windows.Forms.ComboBox cb_fun_act;
        private System.Windows.Forms.Button b_save;
        private System.Windows.Forms.GroupBox groupBox5;
        private System.Windows.Forms.RadioButton rb_solve;
        private System.Windows.Forms.RadioButton rb_net;
        private System.Windows.Forms.Panel p_solve;
        private System.Windows.Forms.RichTextBox rtf_solve;
        private System.Windows.Forms.Button b_open;
        private System.Windows.Forms.Button b_clear;
        private System.Windows.Forms.LinkLabel linkLabel1;
        private System.Windows.Forms.OpenFileDialog ofd_sourseFile;
        private System.Windows.Forms.Label l_filename;
        private System.Windows.Forms.GroupBox groupBox2;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.NumericUpDown nud_q_learn;
        private System.Windows.Forms.Label label25;
        private System.Windows.Forms.TextBox tb_name;
    }
}