% Code to start part HW8 from.
function needleman_wunsch_algo

% A trial sequence with all amino acids:
% ARNDCQEGHILKMFPSTWYV

% The trial sequences in the HW description:
% ARNDCQE
% ARNCDE

% Ask the user for protein sequence 1.
prompt = 'Enter the first amino acid sequence of length 1 to 25: ';
seq1 = input(prompt, 's');
len1 = length(seq1);
fprintf(1, 'Seq 1 %s has %d residues\n', seq1, len1);

% Ask the user for protein sequence 2.
prompt = 'Enter the second amino acid sequence of length 1 to 25: ';
seq2 = input(prompt, 's');
len2 = length(seq2);
fprintf(1, 'Seq 2 %s has %d residues\n', seq2, len2);

% Initialize a 25x25 matrix of zeros.
A = zeros(25,25);
B = zeros(25,25);


% The supplied initMatrix must be modified!

B = initMatrix(B, seq1, seq2);
A = calculateA(A, B, seq1, seq2);
D = A;
matrixD = calculateD(D, seq1, seq2);

% Print the matrix A.
disp('Matrix A');
printMatrix(A, seq1, seq2);

%Print the matrix B
disp('Matrix B');
printMatrix(A, seq1, seq2);

%Print the matrix D
disp('Matrix D')
printD(matrixD, seq1, seq2);

%Trace-back using D matrix
match(matrixD, seq1, seq2);

end

% Print out a matrix
function mat = initMatrix(mat, seq1, seq2)
   numRows = length(seq1) + 1;
   numCols = length(seq2) + 1;
   
   mat(1,1) = 0; % Init the value at (1,1) to 0.
   d=-2;
   
   % Initialize the first column.
   for i=2:numRows 
        mat(i,1) = i*d+2;
   end
   
   % Initialize the first row.
   for c=2:numCols 
        mat(1,c) = c*d+2;
   end
   
   % Init the rest of the matrix.
   for r=2:numRows 
        for c=2:numCols 
            mat(r,c) = scorePair(seq1(r-1),seq2(c-1));
        end
   end
end

%Calculate A
function matA = calculateA(matA, matB, seq1, seq2)
    numRows = length(seq1) + 1;
    numCols = length(seq2) + 1;
    
   
   matA(1,1) = 0; % Init the value at (1,1) to 0.
   gap=-2;
   
   % Initialize the first column.
   for i=2:numRows 
        matA(i,1) = i*gap+2;
   end
   
   % Initialize the first row.
   for c=2:numCols 
        matA(1,c) = c*gap+2;
   end
   
   for r=2:numRows
       for c=2:numCols
           x= matA(r-1,c)+gap;
           y= matA(r-1,c-1) + matB(r,c);
           z= matA(r,c-1)+gap;
           
           matrix = [x,y,z];
           matA(r,c) = max(matrix); 
       end
   end
end

% Print out a matrix
function printMatrix(mat, seq1, seq2)
   len1 = length(seq1);
   len2 = length(seq2);
   
   % Print the first row of sequence letters
   fprintf(1,'\t\t');
   for r=1:len2
    fprintf(1,'%2c\t', seq2(r));
   end
   fprintf(1,'\n');
   
   % Print the first row of scores.
   fprintf(1,'\t');
   for c=1:len2+1 
       fprintf('%2d\t', mat(1,c));
   end
   fprintf(1,'\n');
   
   % Print the rest of the scores.
   for r=2:len1+1
        % Print the first letter of vertical sequence
        fprintf(1, '%c\t', seq1(r-1));
        for c=1:len2+1
            fprintf(1,'%2d\t',mat(r,c));
        end
        fprintf(1,'\n');
   end
end


%Calculate D
function matD = calculateD(matrixD, seq1, seq2)
   numRows = length(seq1) + 1;
   numCols = length(seq2) + 1;
   
   matD(1,1) = 0; % Init the value at (1,1) to 0.
   
   
   % Initialize the first column.
   for i=2:numRows 
        matD(i,1) = 'V';
   end
   
   % Initialize the first row.
   for c=2:numCols 
        matD(1,c) = 'H';
   end
   disp(matD);
   coordinate = matD(1,2);
   fprintf('Coordinate is: %d\n', coordinate);
   for r=2:numRows
       fprintf('Enering Row %d:\n', r);
       for c=2:numCols
          
           x = matrixD(r-1,c); %vertical
           y = matrixD(r-1,c-1); %diagonal
           z = matrixD(r,c-1); %horizontal
           
           matrix = [x,y,z];
           matD(r,c) = max(matrix); 
           
           %Convert to D, H or V
            if y>=x && y>=z
               matD(r,c) = 'D';
            elseif x>z
               matD(r,c) = 'V';
            else
               matD(r,c) = 'H';   
           end
       end
   end
end

% Print out matrix D
function printD(mat, seq1, seq2)
   len1 = length(seq1);
   len2 = length(seq2);
   
   % Print the first row of sequence letters
   fprintf(1,'\t\t');
   for r=1:len2
    fprintf(1,'%2c\t', seq2(r));
   end
   fprintf(1,'\n');
   
   % Print the first row of scores.
   fprintf(1,'\t');
   for c=1:len2+1 
       fprintf('%2c\t', mat(1,c));
   end
   fprintf(1,'\n');
   
   % Print the rest of the scores.
   for r=2:len1+1
        % Print the first letter of vertical sequence
        fprintf(1, '%c\t', seq1(r-1));
        for c=1:len2+1
            fprintf(1,'%2c\t',mat(r,c));
        end
        fprintf(1,'\n');
   end
end

% Find the trace-back alignment of seq1 and seq2
function match(matD, seq1, seq2)
    row = length(seq1) + 1; % Determines number of rows to be looped
    col = length(seq2) + 1; % Determines number of columns to be looped
    
    first_sequence = seq1; % Holds final output of first sequence alignment
    second_sequence = ''; % Holds final ouput of second sequence alignment
    beenVertical = false; % Records if trace-back current position by going vertical
    beenHorizontal = false; % Records if trace-back current position by going horizontal
    
    % Trace through Matrix D backwards (bottom-right to top-left)
    while row>=2 && col>=2

        vertical = matD(row-1,col); % vertical 
        diagonal = matD(row-1,col-1); % diagonal 
        horizontal = matD(row,col-1); % horizontal 
        
        % Checks if diagonal value is the smallest
        if ((diagonal <= horizontal) && (diagonal <= vertical))
        
            % Add to alignment of second sequence
            second_sequence = strcat(second_sequence, seq2(col-1)); 
            % Checks for a mismatch
            if beenVertical == true || beenHorizontal == true
                second_sequence = strcat(second_sequence, '-');
            end 
            % Sets variables back to desired path
            beenVertical = false;
            beenHorizontal = false;
            
            % Moves diagonally in Matrix D
            row = row - 1;
            col = col - 1;
            
        elseif (vertical <= horizontal)         
            % Set variable back to desired path
            beenVertical = true;
            
            % Moves vertically in the Matrix D
            row = row - 1;
       
        else   
            % Set variable back to desired path
            beenHorizontal = true;
           
            % Moves horizontally in the Matrix D
            col = col - 1;
       
        end
        
    end
    
    % Flips the second alignment sequence
    second_sequence = fliplr(second_sequence);
    fprintf('\n');
    % Prints out the alignement of the sequences as calculated
    disp(first_sequence);
    for i=1:length(second_sequence)     
        if first_sequence(i) == second_sequence(i)
            fprintf('|');
        else
            fprintf(' ');
        end
    end
    fprintf('\n');
    disp(second_sequence);
    
end

function value = scorePair(char1, char2)
  % Order: ARNDCQEGHILKMFPSTWYV
  i1 = aaToInt(char1);
  i2 = aaToInt(char2);
  A = [
      4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0;
     -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3;
     -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3;
     -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3;
      0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1;
     -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2;
     -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2;
      0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3;
     -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3;
     -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3;
     -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1;
     -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2;
     -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1;
     -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1;
     -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2;
      1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2;
      0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0;
     -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3;
     -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1;
      0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4];
    value = A(i1, i2);
end

function value = aaToInt(aa)
    switch aa
        case 'A'
            value = 1;
        case 'R' 
            value = 2;
        case 'N'
            value = 3;
        case 'D'
            value = 4;
        case 'C'
            value = 5;
        case 'Q'
            value = 6;
        case 'E'
            value = 7;
        case 'G'
            value = 8;
        case 'H'
            value = 9;
        case 'I'
            value = 10;
        case 'L'
            value = 11;
        case 'K'
            value = 12;
        case 'M'
            value = 13;
        case 'F'
            value = 14;
        case 'P'
            value = 15;
        case 'S'
            value = 16;
        case 'T'
            value = 17;
        case 'W' 
            value = 18;
        case 'Y'
            value = 19;
        case 'V'
            value = 20;
        otherwise
            value = 1;
    end
end