% Reverse output file!
function reverseFile(fileNom)
import t6.*;

%fileNom = 'fwd_1';
theFile = OutputFile(fileNom);
%assert(theFile.numFramesAvailable() == numel(theFile.timesteps()));

timestepBytes = theFile.FrameSize * theFile.BytesPerValue;

maxBufferValues = 1e6;  % number of vals to read in at once
bufferTimesteps = max(1, floor(maxBufferValues / theFile.FrameSize));
buffer = zeros(theFile.FrameSize, bufferTimesteps);

numTimesteps = theFile.numFramesAvailable();
nextFrame = 1;
lastOutFrame = numTimesteps;

fid_in = fopen(theFile.FileName, 'r');
fid_out = fopen([theFile.FileName, '.rev'], 'w');

fileSize = timestepBytes*numTimesteps;
fwrite(fid_out, 0, 'uchar', fileSize-1); % allocate the file a silly way

%%

while nextFrame <= numTimesteps
    framesRemaining = numTimesteps - nextFrame + 1;
    
    framesToRead = min(framesRemaining, bufferTimesteps);
    valuesToRead = framesToRead*theFile.FrameSize;
    
    % Read the frames
    ff = fread(fid_in, valuesToRead, theFile.Precision);
    buffer(1:valuesToRead) = ff(:);
    
    % Reverse the buffer
    buffer = fliplr(buffer(:,1:framesToRead));
    
    % Move the file pointer
    firstOutFrame = lastOutFrame - framesToRead + 1;
    
    BEGINNING_OF_FILE = -1;
    outPos0 = ftell(fid_out);
    status = fseek(fid_out, (firstOutFrame-1)*timestepBytes, BEGINNING_OF_FILE);
    if status == -1
        fprintf('fseek: %s\n', ferror(fid_out));
    end
    outPos = ftell(fid_out);
    
    % Write the buffer
    count = fwrite(fid_out, buffer(1:valuesToRead), theFile.Precision);
    assert(count == valuesToRead);
    
    %fprintf('Out frames %i to %i\n', firstOutFrame, lastOutFrame);
    %fprintf('\tpos is %i\n', outPos);
    %fprintf('Wrote %i values in %i to %i\n', count, outPos0, outPos);
    
    nextFrame = nextFrame + framesToRead;
    lastOutFrame = firstOutFrame-1;
end


fclose(fid_in);
fclose(fid_out);

[pathstr, nom, extension, vers] = fileparts(theFile.SpecFileName);
success = copyfile(theFile.SpecFileName, [nom, '.rev.txt']);

%%

