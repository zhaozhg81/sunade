/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.general;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.nio.ByteBuffer;
import java.util.Date;

/*
 * Based on QRBG.cpp and QRBG.h 
 * 		Designed, written and (c) by Radomir Stevanovic, Jan/2007.
 * 		Developed in Rudjer Boskovic Institute, Zagreb, Croatia.
 * 		Last revision: 2007-Apr-04
 * 
 * Java implementation by Jeroen Rakhorst,The Netherlands, July/2007
 *
 * 
 * You can initialize this class with just a username and a password, or include the server and port to connect to (see constructors)
 * After that, use the readRandomBytes function to read an [length] amount of bytes into the second parameter, the output stream
 * 
 * if you put the static boolean doLog to true, you will see what happens while connecting etc. 
 * Might be useful when having problems while using this class
 * 
 */
public class QRBG {

	public static boolean doLog = false;

	public static enum OperationCodes {
		GET_DATA_AUTH_PLAIN, GET_DATA_AUTH_CERT
	};

	public static enum ServerResponseCodes {
		OK, // everything is ok (user found, quota not exceeded), sending data
		SERVER_STOPPING, // server is stopping (or at least it's shutting down this connection!)
		SERVER_ERROR, // internal server error
		UNKNOWN_OP, // client requested unknown/unsupported operation
		ILL_FORMED_REQUEST, // client sent an ill-formed request packet
		TIMEOUT, // timeout while receiving the request from client
		AUTH_FAILED, // user could not be authenticated - see enum RefusalReasonCodes
		QUOTA_EXCEEDED
		// user quota is (or would be exceeded) - see enum RefusalReasonCodes
	};

	public static enum RefusalReasonCodes {
		NONE(0),

		//
		// bytes per time period quotas
		//
		BYTE_QUOTA_EXCEEDED_FOR_SESSION(0x10), // requested to much data in one session
		BYTE_QUOTA_WOULD_EXCEED_FOR_EON(0x11), // by serving this request, eon quota would be exceeded (request less data)
		BYTE_QUOTA_EXCEEDED_FOR_EON(0x12), // eon quota is already exceeded
		BYTE_QUOTA_WOULD_EXCEED_FOR_YEAR(0x13), // by serving this request, yearly quota would be exceeded (request less data)
		BYTE_QUOTA_EXCEEDED_FOR_YEAR(0x14), // yearly quota is already exceeded
		BYTE_QUOTA_WOULD_EXCEED_FOR_MONTH(0x15), // by serving this request, monthly quota would be exceeded (request less data)
		BYTE_QUOTA_EXCEEDED_FOR_MONTH(0x16), // monthly quota is already exceeded
		BYTE_QUOTA_WOULD_EXCEED_FOR_DAY(0x17), // by serving this request, daily quota would be exceeded (request less data)
		BYTE_QUOTA_EXCEEDED_FOR_DAY(0x18), // daily quota is already exceeded

		//
		// concurrent connections quota AND connection count per time period quotas
		//
		CONCURRENT_CONNECTIONS_QUOTA_EXCEEDED(0x20), // maximum number of allowed parallel requests for authenticated user is already being served (wait and try again)
		CC_QUOTA_EXCEEDED_PER_MINUTE(0x21), // user connections-per-minute limit exceeded (wait and try again)
		CC_QUOTA_EXCEEDED_PER_HOUR(0x22), // user connections-per-hour limit exceeded (wait and try again)
		CC_QUOTA_EXCEEDED_PER_DAY(0x23), // user connections-per-day limit exceeded (wait and try again)
		CC_QUOTA_EXCEEDED_PER_MONTH(0x24), // user connections-per-month limit exceeded (wait and try again)
		CC_QUOTA_EXCEEDED_PER_YEAR(0x25), // user connections-per-year limit exceeded (wait and try again)
		CC_QUOTA_EXCEEDED_PER_EON(0x26) // user connections-per-eon limit exceeded (that's all folks!)
		;
		private final int value;

		RefusalReasonCodes(int value) {
			this.value = value;
		}

		public int getValue() {
			return value;
		}

		public static RefusalReasonCodes getRefusalReasonForValue(int value) {
			for (RefusalReasonCodes r : values()) {
				if (r.getValue() == value) {
					return r;
				}
			}
			return null;
		}
	};

	private ByteBuffer byteBuffer;

	private InetSocketAddress server;

	public QRBG(String username, String password) {
		this("random.irb.hr", 1227, username, password);
	}

	public QRBG(String host, int port, String username, String password) {
		this.server = new InetSocketAddress(host, port);
		int length = 1 + username.length() + 1 + password.length() + 4;
		initializeByteBuffer(OperationCodes.GET_DATA_AUTH_PLAIN, username, password, length);
	}

	/*
	 * This function is a small preparation for the possible use of secure communication 
	 */
	private void initializeByteBuffer(OperationCodes operationCode, String username, String password, int length) {
		switch (operationCode) {
		case GET_DATA_AUTH_PLAIN:
			byteBuffer = ByteBuffer.allocate(length + 3);
			byteBuffer.put((byte) operationCode.ordinal());
			byteBuffer.putShort(new Integer(length).shortValue());
			byteBuffer.put((byte) username.length());
			byteBuffer.put(username.getBytes());
			byteBuffer.put((byte) password.length());
			byteBuffer.put(password.getBytes());
			break;
		default:
			throw new UnsupportedOperationException("Cannot handle operation of type " + operationCode);
		}
	}

	/*
	 * connect with server and if the server returns random bytes, write them to the outputstream
	 * @returns the amount of bytes actually written to the outputstream
	 */
	public int readRandomBytes(int length, OutputStream output) throws QRBG_Exception, IOException {
		byteBuffer.putInt(byteBuffer.capacity() - 4, length);
		log("opening socket");
		int actuallyWritten = 0;
		Socket socket = new Socket();
		try {
			socket.connect(this.server);
			log("connected to server " + this.server);
			OutputStream serverOutputStream = socket.getOutputStream();
			serverOutputStream.write(byteBuffer.array());
			log(byteBuffer.capacity() + " bytes written to server");
			serverOutputStream.flush();
			InputStream serverInputStream = socket.getInputStream();
			int result;
			log("reading from server");
			boolean readServerStatus = false;
			boolean readReplySize = false;
			int resultingSize = 0;
			while ((result = serverInputStream.read()) != -1) {
				if (!readServerStatus) {
					processResponseCodes(serverInputStream, result);
					readServerStatus = true;
				} else if (!readReplySize) {
					resultingSize = processReplySize(serverInputStream, result);
					log("server answer to amount of returning random bytes : " + resultingSize);
					readReplySize = true;
				} else {
					output.write(result);
					actuallyWritten++;
				}
			}
			log("Finished reading " + actuallyWritten + " bytes");
		} finally {
			if (socket.isConnected()) {
				try {
					socket.shutdownOutput();
					socket.shutdownInput();
					socket.close();
				} catch (IOException e) {
					log(e);
				}
			}
		}
		return actuallyWritten;
	}

	private void log(String message) {
		if (doLog) {
			System.out.println(new Date() + " - " + message);
		}
	}

	private void log(Exception e) {
		if (doLog) {
			e.printStackTrace();
		}
	}

	private int processReplySize(InputStream inputStream, int byteValue) throws IOException {
		ByteBuffer b = ByteBuffer.allocate(4);
		b.put((byte) byteValue);
		b.put(new Integer(inputStream.read()).byteValue());
		b.put(new Integer(inputStream.read()).byteValue());
		b.put(new Integer(inputStream.read()).byteValue());
		return b.getInt(0);
	}

	private void processResponseCodes(InputStream inputStream, int byteValue) throws QRBG_Exception, IOException {
		int result;
		ServerResponseCodes responseCodes = ServerResponseCodes.values()[byteValue];
		result = inputStream.read();
		byteValue = new Integer(result).byteValue();
		if (responseCodes != ServerResponseCodes.OK) {
			throw new QRBG_Exception(responseCodes, RefusalReasonCodes.getRefusalReasonForValue(byteValue));
		}
	}

	public class QRBG_Exception extends Exception {
		private static final long serialVersionUID = 5515425944581127871L;

		private final ServerResponseCodes serverResponseCodes;

		private final RefusalReasonCodes refusalReasonCodes;

		public QRBG_Exception(ServerResponseCodes serverResponseCodes, RefusalReasonCodes refusalReasonCodes) {
			this.serverResponseCodes = serverResponseCodes;
			this.refusalReasonCodes = refusalReasonCodes;

		}

		public RefusalReasonCodes getRefusalReasonCodes() {
			return refusalReasonCodes;
		}

		public ServerResponseCodes getServerResponseCodes() {
			return serverResponseCodes;
		}

		@Override
		public String getMessage() {
			return "Error processing request : "+serverResponseCodes.toString()+", reason : "+refusalReasonCodes.toString();
		}
	}
	
	public static void main(String[] args) throws IOException
	{
            /*		if (args.length<3)
		{
			System.out.println("Usage : QRBG [username] [password] [amountOfBytesToRead]");
			return;
		}
            */  
		QRBG qrbg=new QRBG("zhaozhg","z19810329");
		Integer numBytes=Integer.parseInt("100");
		java.io.ByteArrayOutputStream baos=new java.io.ByteArrayOutputStream(numBytes);
		try {
			System.out.println("resulting bytes : "+qrbg.readRandomBytes(numBytes,baos));
			int lineLength=20;
			for (byte b:baos.toByteArray())
			{
				System.out.print(b);
				System.out.print(" ");
				lineLength--;
				if (lineLength<=0)
				{
					System.out.println();
					lineLength=20;
				}
			}
		} catch (QRBG_Exception e) {
			System.err.println("Error reading bytes from server, server response : "+e.getMessage());
		}
	}
}

