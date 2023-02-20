LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;

ENTITY DIVIDER IS
GENERIC
(
	DIVIDER_WIDTH: INTEGER := 16;
	TUNING_FACTOR: INTEGER := 6
);
PORT
(
	CLK, RST: IN STD_LOGIC;
	Z		: IN STD_LOGIC_VECTOR(DIVIDER_WIDTH-1 DOWNTO 0);
	D		: IN STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0);
	Q		: OUT STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0);
	S		: OUT STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0)
);
END ENTITY DIVIDER;

ARCHITECTURE STRUCTURAL OF DIVIDER IS
-- INTERNAL SIGNALS FOR CELL CONNECTION
SIGNAL Z_TEMP: STD_LOGIC_VECTOR(DIVIDER_WIDTH-1 DOWNTO 0);
SIGNAL D_TEMP: STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0);
TYPE TEMP_REG_SIGNAL IS ARRAY (0 TO (DIVIDER_WIDTH/2)-1) OF STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0);
SIGNAL BOUT_TEMP, DIFF_TEMP, S_TEMP, Q_TEMP: TEMP_REG_SIGNAL; 

-- SIGNAL DECLARATION FOR PIPELINING
TYPE PIPE_REG_SIGNAL IS ARRAY (0 TO (DIVIDER_WIDTH/2)-1) OF STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-1 DOWNTO 0);
SIGNAL Q_PIPE_OUT, S_PIPE_OUT, DIVS_PIPE_OUT: PIPE_REG_SIGNAL;

TYPE DIVD_PIPE_REG_SIGNAL IS ARRAY (0 TO (DIVIDER_WIDTH/2)-1) OF STD_LOGIC_VECTOR((DIVIDER_WIDTH/2)-2 DOWNTO 0);
SIGNAL DIVD_PIPE_OUT: DIVD_PIPE_REG_SIGNAL;

BEGIN
	-- ASSIGN INPUTS
	Z_TEMP <= Z;
	D_TEMP <= D;
	
	ZERO_L0_TUNING: IF (TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) <= 0 GENERATE
		EX_Q_L0_0: ENTITY WORK.EXACT_CELL_Q PORT MAP
			(
				X 	 => Z_TEMP((DIVIDER_WIDTH/2)-1), 
				Y 	 => D_TEMP(0), BIN => '0', 
				D 	 => DIFF_TEMP(0)(0), 
				BOUT => BOUT_TEMP(0)(0)
			);
			
		STAGE_Q_L0_0: FOR J IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			EX_Q_J: ENTITY WORK.EXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(J+((DIVIDER_WIDTH/2)-1)),
				Y    => D_TEMP(J),
				BIN  => BOUT_TEMP(0)(J-1),
				D    => DIFF_TEMP(0)(J),
				BOUT => BOUT_TEMP(0)(J)
			);
			END GENERATE STAGE_Q_L0_0;

		Q_TEMP(0)((DIVIDER_WIDTH/2)-1) <= (Z_TEMP(DIVIDER_WIDTH-1) OR NOT(BOUT_TEMP(0)((DIVIDER_WIDTH/2)-1)));
		Q_TEMP(0)((DIVIDER_WIDTH/2)-2 downto 0) <= (others => '0');

		STAGE_R_L0_0: FOR R IN 0 TO DIVIDER_WIDTH/2-1 GENERATE
			EX_R_R: ENTITY WORK.EXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(R+(DIVIDER_WIDTH/2-1)),
				D  => DIFF_TEMP(0)(R),
				R  => S_TEMP(0)(R)
			);
			END GENERATE STAGE_R_L0_0;
	END GENERATE ZERO_L0_TUNING;
		
	ZERO_L1_TUNING: IF (TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) = 1 GENERATE
		IEX_Q_L1_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
			(
				X 	 => Z_TEMP((DIVIDER_WIDTH/2)-1), 
				Y 	 => D_TEMP(0), 
				BIN => '0', 
				D 	 => DIFF_TEMP(0)(0), 
				BOUT => BOUT_TEMP(0)(0)
			);
				
		STAGE_Q_L1_0: FOR J IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			EX_Q_J: ENTITY WORK.EXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(J+((DIVIDER_WIDTH/2)-1)),
				Y    => D_TEMP(J),
				BIN  => BOUT_TEMP(0)(J-1),
				D    => DIFF_TEMP(0)(J),
				BOUT => BOUT_TEMP(0)(J)
			);
			END GENERATE STAGE_Q_L1_0;

		Q_TEMP(0)((DIVIDER_WIDTH/2)-1) <= (Z_TEMP(DIVIDER_WIDTH-1) OR NOT(BOUT_TEMP(0)((DIVIDER_WIDTH/2)-1)));
		Q_TEMP(0)((DIVIDER_WIDTH/2)-2 downto 0) <= (others => '0');
		
		IEX_R_0: ENTITY WORK.INEXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1), 
				X  => Z_TEMP(DIVIDER_WIDTH/2-1), 
				D  => DIFF_TEMP(0)(0), 
				R  => S_TEMP(0)(0)
			);

		STAGE_R_L1_0: FOR R IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			EX_R_R: ENTITY WORK.EXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(R+(DIVIDER_WIDTH/2-1)),
				D  => DIFF_TEMP(0)(R),
				R  => S_TEMP(0)(R)
			);
			END GENERATE STAGE_R_L1_0;
	END GENERATE ZERO_L1_TUNING;
	
	ZERO_L2_TUNING: IF ((TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) > 1 and (TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) < DIVIDER_WIDTH/2-1) GENERATE
		IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
			(
				X 	 => Z_TEMP((DIVIDER_WIDTH/2)-1), 
				Y 	 => D_TEMP(0), BIN => '0', 
				D 	 => DIFF_TEMP(0)(0), 
				BOUT => BOUT_TEMP(0)(0)
			);
		
		STAGE_IEQ_L2_0: FOR J IN 1 TO (TUNING_FACTOR - DIVIDER_WIDTH/2) GENERATE
			IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(J+((DIVIDER_WIDTH/2)-1)),
				Y    => D_TEMP(J),
				BIN  => BOUT_TEMP(0)(J-1),
				D    => DIFF_TEMP(0)(J),
				BOUT => BOUT_TEMP(0)(J)
			);
			END GENERATE STAGE_IEQ_L2_0;
			
		STAGE_EQ_L2_0: FOR K IN (TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) TO (DIVIDER_WIDTH/2-1) GENERATE
			IEX_Q_K: ENTITY WORK.EXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(K+((DIVIDER_WIDTH/2)-1)),
				Y    => D_TEMP(K),
				BIN  => BOUT_TEMP(0)(K-1),
				D    => DIFF_TEMP(0)(K),
				BOUT => BOUT_TEMP(0)(K)
			);
			END GENERATE STAGE_EQ_L2_0;
				
		Q_TEMP(0)((DIVIDER_WIDTH/2)-1) <= (Z_TEMP(DIVIDER_WIDTH-1) OR NOT(BOUT_TEMP(0)((DIVIDER_WIDTH/2)-1)));
		Q_TEMP(0)((DIVIDER_WIDTH/2)-2 downto 0) <= (others => '0');
			
		STAGE_IER_0: FOR R IN 0 TO (TUNING_FACTOR - DIVIDER_WIDTH/2) GENERATE
			IEX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(R+(DIVIDER_WIDTH/2-1)),
				D  => DIFF_TEMP(0)(R),
				R  => S_TEMP(0)(R)
			);
			END GENERATE STAGE_IER_0;
			
		STAGE_ER_0: FOR S IN (TUNING_FACTOR + 1 - DIVIDER_WIDTH/2) TO DIVIDER_WIDTH/2-1 GENERATE
			EX_R_S: ENTITY WORK.EXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(S+(DIVIDER_WIDTH/2-1)),
				D  => DIFF_TEMP(0)(S),
				R  => S_TEMP(0)(S)
			);
			END GENERATE STAGE_ER_0;
	END GENERATE ZERO_L2_TUNING;
	
	ZERO_L3_TUNING: IF ((TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) = (DIVIDER_WIDTH/2-1)) GENERATE
		IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
			(
				X 	 => Z_TEMP((DIVIDER_WIDTH/2)-1), 
				Y 	 => D_TEMP(0), BIN => '0', 
				D 	 => DIFF_TEMP(0)(0), 
				BOUT => BOUT_TEMP(0)(0)
			);
		
		STAGE_IEQ_L3_0: FOR J IN 1 TO (TUNING_FACTOR - DIVIDER_WIDTH/2) GENERATE
			IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(J+((DIVIDER_WIDTH/2)-1)),
				Y    => D_TEMP(J),
				BIN  => BOUT_TEMP(0)(J-1),
				D    => DIFF_TEMP(0)(J),
				BOUT => BOUT_TEMP(0)(J)
			);
			END GENERATE STAGE_IEQ_L3_0;
			
			EX_Q_7: ENTITY WORK.EXACT_CELL_Q PORT MAP
			(
				X    => Z_TEMP(TUNING_FACTOR),
				Y    => D_TEMP(TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)),
				BIN  => BOUT_TEMP(0)(TUNING_FACTOR - (DIVIDER_WIDTH/2)),
				D    => DIFF_TEMP(0)(TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)),
				BOUT => BOUT_TEMP(0)(TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2))
			);
				
		Q_TEMP(0)((DIVIDER_WIDTH/2)-1) <= (Z_TEMP(DIVIDER_WIDTH-1) OR NOT(BOUT_TEMP(0)((DIVIDER_WIDTH/2)-1)));
		Q_TEMP(0)((DIVIDER_WIDTH/2)-2 downto 0) <= (others => '0');
			
		STAGE_IER_0: FOR R IN 0 TO (TUNING_FACTOR - DIVIDER_WIDTH/2) GENERATE
			IEX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(R+(DIVIDER_WIDTH/2-1)),
				D  => DIFF_TEMP(0)(R),
				R  => S_TEMP(0)(R)
			);
			END GENERATE STAGE_IER_0;
			
			EX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
			(
				QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
				X  => Z_TEMP(TUNING_FACTOR),
				D  => DIFF_TEMP(0)(TUNING_FACTOR + 1 - DIVIDER_WIDTH/2),
				R  => S_TEMP(0)(TUNING_FACTOR + 1 - DIVIDER_WIDTH/2)
			);
	END GENERATE ZERO_L3_TUNING;

	ZERO_L4_TUNING: IF (TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) >= DIVIDER_WIDTH/2 GENERATE
			IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
					X 	 => Z_TEMP((DIVIDER_WIDTH/2)-1), 
					Y 	 => D_TEMP(0), BIN => '0', 
					D 	 => DIFF_TEMP(0)(0), 
					BOUT => BOUT_TEMP(0)(0)
				);
			
			STAGE_IEQ_0: FOR J IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
				IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
					X    => Z_TEMP(J+((DIVIDER_WIDTH/2)-1)),
					Y    => D_TEMP(J),
					BIN  => BOUT_TEMP(0)(J-1),
					D    => DIFF_TEMP(0)(J),
					BOUT => BOUT_TEMP(0)(J)
				);
				END GENERATE STAGE_IEQ_0;
					
			Q_TEMP(0)((DIVIDER_WIDTH/2)-1) <= (Z_TEMP(DIVIDER_WIDTH-1) OR NOT(BOUT_TEMP(0)((DIVIDER_WIDTH/2)-1)));
			Q_TEMP(0)((DIVIDER_WIDTH/2)-2 downto 0) <= (others => '0');
				
			STAGE_IER_0: FOR R IN 0 TO DIVIDER_WIDTH/2-1 GENERATE
				IEX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
					QS => Q_TEMP(0)(DIVIDER_WIDTH/2-1),
					X  => Z_TEMP(R+(DIVIDER_WIDTH/2-1)),
					D  => DIFF_TEMP(0)(R),
					R  => S_TEMP(0)(R)
				);
				END GENERATE STAGE_IER_0;	
	END GENERATE ZERO_L4_TUNING;
	
	-- Stage 0 Pipeline Registers
	REG_Q_S_0: ENTITY WORK.REG 
		GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
		PORT MAP
		(
			CLK 	=> CLK,
		    RST 	=> RST,
		    INPUT   => Q_TEMP(0),
		    OUTPUT 	=> Q_PIPE_OUT(0)
		);
	
	REG_S_S_0: ENTITY WORK.REG 
		GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
		PORT MAP
		(
			CLK 	=> CLK,
		    RST 	=> RST,
		    INPUT 	=> S_TEMP(0),
		    OUTPUT 	=> S_PIPE_OUT(0)
		);
	
	REG_DS_S_0: ENTITY WORK.REG 
		GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
		PORT MAP
		(
			CLK 	=> CLK,
		    RST 	=> RST,
		    INPUT  	=> D_TEMP,
		    OUTPUT 	=> DIVS_PIPE_OUT(0)
		);
		
	REG_DD_S_0: ENTITY WORK.REG 
		GENERIC MAP (INPUT_WIDTH => (DIVIDER_WIDTH/2 - 1))
		PORT MAP
		(
			CLK 	=> CLK,
		    RST 	=> RST,
		    INPUT  	=> Z_TEMP((DIVIDER_WIDTH/2 - 2) DOWNTO 0),
		    OUTPUT 	=> DIVD_PIPE_OUT(0)
		);
	
	ROW_GEN: FOR I IN 1 TO (DIVIDER_WIDTH/2)-1 GENERATE
		QNT_L1: IF I = 1 GENERATE
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-1) <= Q_PIPE_OUT(I-1)((DIVIDER_WIDTH/2)-1);
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-2) <= (S_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1) OR NOT(BOUT_TEMP(I)(DIVIDER_WIDTH/2-1)));
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-3 downto 0) <= (others => '0');
		END GENERATE;
		
		QNT_L2: IF I = (DIVIDER_WIDTH/2)-2 GENERATE
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-1 DOWNTO 2) <= Q_PIPE_OUT(I-1)((DIVIDER_WIDTH/2)-1 DOWNTO 2);
		   Q_TEMP(I)(1) <= (S_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1) OR NOT(BOUT_TEMP(I)(DIVIDER_WIDTH/2-1)));
		   Q_TEMP(I)(0) <= '0';
		END GENERATE;
		
		QNT_L3: IF I = (DIVIDER_WIDTH/2)-1 GENERATE
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-1 DOWNTO 1) <= Q_PIPE_OUT(I-1)((DIVIDER_WIDTH/2)-1 DOWNTO 1);
		   Q_TEMP(I)(0) <= (S_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1) OR NOT(BOUT_TEMP(I)(DIVIDER_WIDTH/2-1)));
		END GENERATE;
		
		QNT_L4: IF (I > 1 AND I < (DIVIDER_WIDTH/2)-2) GENERATE
		   Q_TEMP(I)((DIVIDER_WIDTH/2)-1 DOWNTO (DIVIDER_WIDTH/2) - I) <= Q_PIPE_OUT(I-1)((DIVIDER_WIDTH/2)-1 DOWNTO (DIVIDER_WIDTH/2) - I);
		   Q_TEMP(I)((DIVIDER_WIDTH/2)- 1 - I) <= (S_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1) OR NOT(BOUT_TEMP(I)(DIVIDER_WIDTH/2-1)));
		   Q_TEMP(I)((DIVIDER_WIDTH/2)- 2 - I downto 0) <= (others => '0');
		END GENERATE;
		
		-- OTHER Stage Pipeline Registers
		REG_Q_S_I: ENTITY WORK.REG 
			GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
			PORT MAP
			(
				CLK 	=> CLK,
				RST 	=> RST,
				INPUT   => Q_TEMP(I),
				OUTPUT 	=> Q_PIPE_OUT(I)
			);
		
		REG_S_S_I: ENTITY WORK.REG 
			GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
			PORT MAP
			(
				CLK 	=> CLK,
				RST 	=> RST,
				INPUT 	=> S_TEMP(I),
				OUTPUT 	=> S_PIPE_OUT(I)
			);
		
		REG_DS_S_I: ENTITY WORK.REG 
			GENERIC MAP (INPUT_WIDTH => DIVIDER_WIDTH/2)
			PORT MAP
			(
				CLK 	=> CLK,
				RST 	=> RST,
				INPUT  	=> DIVS_PIPE_OUT(I-1),
				OUTPUT 	=> DIVS_PIPE_OUT(I)
			);
			
		REG_DD_S_I: ENTITY WORK.REG 
			GENERIC MAP (INPUT_WIDTH => (DIVIDER_WIDTH/2 - 1))
			PORT MAP
			(
				CLK 	=> CLK,
				RST 	=> RST,
				INPUT  	=> DIVD_PIPE_OUT(I-1),
				OUTPUT 	=> DIVD_PIPE_OUT(I)
			);
		
	    OTHER_L0_TUNING: IF (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) <= 0 GENERATE
			EX_Q_0: ENTITY WORK.EXACT_CELL_Q PORT MAP
				(
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					Y => DIVS_PIPE_OUT(I-1)(0), 
					BIN => '0', 
					D => DIFF_TEMP(I)(0), 
					BOUT => BOUT_TEMP(I)(0)
				);
		
			STAGE_Q_I: FOR J IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			  EX_Q_J: ENTITY WORK.EXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(J-1),
				  Y    => DIVS_PIPE_OUT(I-1)(J),
				  BIN  => BOUT_TEMP(I)(J-1),
				  D    => DIFF_TEMP(I)(J),
				  BOUT => BOUT_TEMP(I)(J)
				);
			END GENERATE STAGE_Q_I;
			
			
			EX_R_0: ENTITY WORK.EXACT_CELL_R PORT MAP
				(
					QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I), 
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					D => DIFF_TEMP(I)(0), 
					R => S_TEMP(I)(0)
				);
		  
			STAGE_R_I: FOR R IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			  EX_R_R: ENTITY WORK.EXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(R-1),
				  D  => DIFF_TEMP(I)(R),
				  R  => S_TEMP(I)(R)
				);
			END GENERATE STAGE_R_I;	
	    END GENERATE OTHER_L0_TUNING;
		
		OTHER_L1_TUNING: IF ((I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) > 0 and (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) < (DIVIDER_WIDTH/2-1)) GENERATE
			IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					Y => DIVS_PIPE_OUT(I-1)(0), 
					BIN => '0', 
					D => DIFF_TEMP(I)(0), 
					BOUT => BOUT_TEMP(I)(0)
				);
		
			STAGE_IEQ_I: FOR J IN 1 TO (I + TUNING_FACTOR - (DIVIDER_WIDTH/2)) GENERATE
			  IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(J-1),
				  Y    => DIVS_PIPE_OUT(I-1)(J),
				  BIN  => BOUT_TEMP(I)(J-1),
				  D    => DIFF_TEMP(I)(J),
				  BOUT => BOUT_TEMP(I)(J)
				);
			END GENERATE STAGE_IEQ_I;
			
			STAGE_EQ_I: FOR K IN (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) TO DIVIDER_WIDTH/2-1 GENERATE
			  EX_Q_K: ENTITY WORK.EXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(K-1),
				  Y    => DIVS_PIPE_OUT(I-1)(K),
				  BIN  => BOUT_TEMP(I)(K-1),
				  D    => DIFF_TEMP(I)(K),
				  BOUT => BOUT_TEMP(I)(K)
				);
			END GENERATE STAGE_EQ_I;
				
			IEX_R_0: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
					QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I), 
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					D => DIFF_TEMP(I)(0), 
					R => S_TEMP(I)(0)
				);
		  
			STAGE_IER_I: FOR R IN 1 TO (I + TUNING_FACTOR - (DIVIDER_WIDTH/2)) GENERATE
			  IEX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(R-1),
				  D  => DIFF_TEMP(I)(R),
				  R  => S_TEMP(I)(R)
				);
			END GENERATE STAGE_IER_I;	
			
			STAGE_ER_I: FOR S IN (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) TO DIVIDER_WIDTH/2-1 GENERATE
			  EX_R_S: ENTITY WORK.EXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(S-1),
				  D  => DIFF_TEMP(I)(S),
				  R  => S_TEMP(I)(S)
				);
			END GENERATE STAGE_ER_I;
	    END GENERATE OTHER_L1_TUNING;
		
		OTHER_L2_TUNING: IF (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2) = (DIVIDER_WIDTH/2-1)) GENERATE
			IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					Y => DIVS_PIPE_OUT(I-1)(0), 
					BIN => '0', 
					D => DIFF_TEMP(I)(0), 
					BOUT => BOUT_TEMP(I)(0)
				);
		
			STAGE_IEQ_I: FOR J IN 1 TO (I + TUNING_FACTOR - (DIVIDER_WIDTH/2)) GENERATE
			  IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(J-1),
				  Y    => DIVS_PIPE_OUT(I-1)(J),
				  BIN  => BOUT_TEMP(I)(J-1),
				  D    => DIFF_TEMP(I)(J),
				  BOUT => BOUT_TEMP(I)(J)
				);
			END GENERATE STAGE_IEQ_I;
			
			  EX_Q_K: ENTITY WORK.EXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(I + TUNING_FACTOR - (DIVIDER_WIDTH/2)),
				  Y    => DIVS_PIPE_OUT(I-1)(I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)),
				  BIN  => BOUT_TEMP(I)(I + TUNING_FACTOR - (DIVIDER_WIDTH/2)),
				  D    => DIFF_TEMP(I)(I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)),
				  BOUT => BOUT_TEMP(I)(I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2))
				);
				
			IEX_R_0: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
					QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I), 
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					D => DIFF_TEMP(I)(0), 
					R => S_TEMP(I)(0)
				);
		  
			STAGE_IER_I: FOR R IN 1 TO (I + TUNING_FACTOR - (DIVIDER_WIDTH/2)) GENERATE
			  IEX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(R-1),
				  D  => DIFF_TEMP(I)(R),
				  R  => S_TEMP(I)(R)
				);
			END GENERATE STAGE_IER_I;	
			
			  EX_R_S: ENTITY WORK.EXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(I + TUNING_FACTOR - (DIVIDER_WIDTH/2)),
				  D  => DIFF_TEMP(I)(I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)),
				  R  => S_TEMP(I)(I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2))
				);
	    END GENERATE OTHER_L2_TUNING;
		
		OTHER_L3_TUNING: IF (I + TUNING_FACTOR + 1 - (DIVIDER_WIDTH/2)) >= DIVIDER_WIDTH/2 GENERATE
			IEX_Q_0: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					Y => DIVS_PIPE_OUT(I-1)(0), 
					BIN => '0', 
					D => DIFF_TEMP(I)(0), 
					BOUT => BOUT_TEMP(I)(0)
				);
		
			STAGE_IEQ_I: FOR J IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			  IEX_Q_J: ENTITY WORK.INEXACT_CELL_Q PORT MAP
				(
				  X    => S_PIPE_OUT(I-1)(J-1),
				  Y    => DIVS_PIPE_OUT(I-1)(J),
				  BIN  => BOUT_TEMP(I)(J-1),
				  D    => DIFF_TEMP(I)(J),
				  BOUT => BOUT_TEMP(I)(J)
				);
			END GENERATE STAGE_IEQ_I;
				
			IEX_R_0: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
					QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I), 
					X => DIVD_PIPE_OUT(I-1)(DIVIDER_WIDTH/2-1-I), 
					D => DIFF_TEMP(I)(0), 
					R => S_TEMP(I)(0)
				);
		  
			STAGE_IER_I: FOR R IN 1 TO DIVIDER_WIDTH/2-1 GENERATE
			  EX_R_R: ENTITY WORK.INEXACT_CELL_R PORT MAP
				(
				  QS => Q_TEMP(I)(DIVIDER_WIDTH/2-1-I),
				  X  => S_PIPE_OUT(I-1)(R-1),
				  D  => DIFF_TEMP(I)(R),
				  R  => S_TEMP(I)(R)
				);
			END GENERATE STAGE_IER_I;	
	    END GENERATE OTHER_L3_TUNING;
	END GENERATE ROW_GEN;
			
	-- ASSIGNING OUTPUTS
	Q <= Q_PIPE_OUT(DIVIDER_WIDTH/2-1);
	S <= S_PIPE_OUT(DIVIDER_WIDTH/2-1);
	
END ARCHITECTURE STRUCTURAL;

		   
		     

