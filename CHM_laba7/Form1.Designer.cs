namespace SLAU
{
    partial class Form1
    {
        /// <summary>
        /// Обязательная переменная конструктора.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Освободить все используемые ресурсы.
        /// </summary>
        /// <param name="disposing">истинно, если управляемый ресурс должен быть удален; иначе ложно.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Код, автоматически созданный конструктором форм Windows

        /// <summary>
        /// Требуемый метод для поддержки конструктора — не изменяйте 
        /// содержимое этого метода с помощью редактора кода.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(Form1));
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.SolveBox = new System.Windows.Forms.ListBox();
            this.groupBox2 = new System.Windows.Forms.GroupBox();
            this.SolveButton = new System.Windows.Forms.Button();
            this.rb_13 = new System.Windows.Forms.RadioButton();
            this.rb_22 = new System.Windows.Forms.RadioButton();
            this.groupBox1.SuspendLayout();
            this.groupBox2.SuspendLayout();
            this.SuspendLayout();
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.SolveBox);
            this.groupBox1.Location = new System.Drawing.Point(13, 13);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(675, 161);
            this.groupBox1.TabIndex = 0;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "Системные сообщения";
            // 
            // SolveBox
            // 
            this.SolveBox.FormattingEnabled = true;
            this.SolveBox.Location = new System.Drawing.Point(7, 20);
            this.SolveBox.Name = "SolveBox";
            this.SolveBox.Size = new System.Drawing.Size(661, 134);
            this.SolveBox.TabIndex = 0;
            // 
            // groupBox2
            // 
            this.groupBox2.Controls.Add(this.rb_22);
            this.groupBox2.Controls.Add(this.rb_13);
            this.groupBox2.Controls.Add(this.SolveButton);
            this.groupBox2.Location = new System.Drawing.Point(694, 13);
            this.groupBox2.Name = "groupBox2";
            this.groupBox2.Size = new System.Drawing.Size(200, 161);
            this.groupBox2.TabIndex = 1;
            this.groupBox2.TabStop = false;
            this.groupBox2.Text = "Панель управления";
            // 
            // SolveButton
            // 
            this.SolveButton.Enabled = false;
            this.SolveButton.Font = new System.Drawing.Font("Microsoft Sans Serif", 9.75F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.SolveButton.Location = new System.Drawing.Point(7, 93);
            this.SolveButton.Name = "SolveButton";
            this.SolveButton.Size = new System.Drawing.Size(187, 61);
            this.SolveButton.TabIndex = 1;
            this.SolveButton.Text = "Решить и сохранить решение";
            this.SolveButton.UseVisualStyleBackColor = true;
            this.SolveButton.Click += new System.EventHandler(this.SolveButton_Click);
            // 
            // rb_13
            // 
            this.rb_13.AutoSize = true;
            this.rb_13.Checked = true;
            this.rb_13.Location = new System.Drawing.Point(27, 33);
            this.rb_13.Name = "rb_13";
            this.rb_13.Size = new System.Drawing.Size(74, 17);
            this.rb_13.TabIndex = 2;
            this.rb_13.TabStop = true;
            this.rb_13.Text = "Номер 13";
            this.rb_13.UseVisualStyleBackColor = true;
            // 
            // rb_22
            // 
            this.rb_22.AutoSize = true;
            this.rb_22.Location = new System.Drawing.Point(27, 56);
            this.rb_22.Name = "rb_22";
            this.rb_22.Size = new System.Drawing.Size(74, 17);
            this.rb_22.TabIndex = 3;
            this.rb_22.Text = "Номер 22";
            this.rb_22.UseVisualStyleBackColor = true;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(897, 183);
            this.Controls.Add(this.groupBox2);
            this.Controls.Add(this.groupBox1);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedDialog;
           // this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.MaximizeBox = false;
            this.MinimizeBox = false;
            this.Name = "Form1";
            this.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
            this.Text = "Решение СЛАУ через LU-разложение";
            this.groupBox1.ResumeLayout(false);
            this.groupBox2.ResumeLayout(false);
            this.groupBox2.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.GroupBox groupBox2;
        public System.Windows.Forms.ListBox SolveBox;
        private System.Windows.Forms.Button SolveButton;
        private System.Windows.Forms.RadioButton rb_22;
        private System.Windows.Forms.RadioButton rb_13;
    }
}

